// Metrics and charts

import { useEffect, useMemo, useState } from "react";
import { MetricDefinitions, MetricHistoryItem, MetricsDashboardProps, MetricsDashboardState, TreeNode, VisibleMetrics } from "../types";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';


// Metric definitions for extensibility
const metricDefinitions: MetricDefinitions = {
    cost: { 
        label: 'Reaction Cost ($)', 
        color: '#EC4899',
        calculate: (nodes: TreeNode[]) => nodes.reduce((sum, node) => sum + (node.cost || 0), 0)
    },
    bandgap: { 
        label: 'Band Gap (eV)', 
        color: '#F59E0B',
        calculate: (nodes: TreeNode[]) => nodes[nodes.length-1]?.bandgap || 0
    },
    density: { 
        label: 'Molecular Density (g/cm³)', 
        color: '#10B981',
        calculate: (nodes: TreeNode[]) => nodes[nodes.length-1]?.density || 0
    },
    /*yield: { 
        label: 'Yield (%)', 
        color: '#10B981',
        calculate: (nodes) => nodes.length > 0 ? (nodes.reduce((sum, node) => sum + (node.yield || Math.random() * 100), 0) / nodes.length) : 0
    },*/
};


// Calculate all metrics
const calculateMetrics = (nodes: TreeNode[]): MetricHistoryItem => {
    const metrics: MetricHistoryItem = { nodeCount: nodes.length, step: 0 };
    Object.keys(metricDefinitions).forEach(key => {
        metrics[key] = metricDefinitions[key].calculate(nodes);
    });
    return metrics;
};

export const useMetricsDashboardState = (): MetricsDashboardState => {
    const [metricsHistory, setMetricsHistory] = useState<MetricHistoryItem[]>([]);
    const [visibleMetrics, setVisibleMetrics] = useState<VisibleMetrics>({
        cost: true,
        bandgap: false,
        density: false,
        yield: false,
    });

    return {metricsHistory, setMetricsHistory, visibleMetrics, setVisibleMetrics};
};

export const MetricsDashboard: React.FC<MetricsDashboardProps> = ({treeNodes, metricsHistory, setMetricsHistory, visibleMetrics, setVisibleMetrics}) => {
    // Memoize the metrics charts to prevent re-render on mouse move
    const metricsCharts = useMemo(() => {
        if (metricsHistory.length === 0) return null;

        return (
            <div className={`grid gap-6 ${Object.values(visibleMetrics).filter(Boolean).length > 1 ? 'grid-cols-2' : 'grid-cols-1'}`}>
            {Object.keys(metricDefinitions).filter(key => visibleMetrics[key as keyof VisibleMetrics]).map(metricKey => {
                const metric = metricDefinitions[metricKey];
                return (
                <div key={metricKey} className="bg-white/5 rounded-xl p-4">
                    <h4 className="text-sm font-semibold text-purple-200 mb-2">{metric.label}</h4>
                    <ResponsiveContainer width="100%" height={200}>
                    <LineChart data={metricsHistory}>
                        <CartesianGrid strokeDasharray="3 3" stroke="#8B5CF6" opacity={0.1} />
                        <XAxis 
                        dataKey="step" 
                        stroke="#A78BFA"
                        tick={{ fontSize: 12 }}
                        />
                        <YAxis 
                        stroke="#A78BFA"
                        tick={{ fontSize: 12 }}
                        />
                        <Tooltip 
                        contentStyle={{ 
                            backgroundColor: '#1E1B4B', 
                            border: '2px solid #8B5CF6',
                            borderRadius: '8px',
                            color: '#E9D5FF'
                        }}
                        // NOTE: Returning `as any` since no other type information worked
                        formatter={(value: number) => value.toFixed(2) as any}
                        />
                        <Line 
                        type="monotone" 
                        dataKey={metricKey}
                        stroke={metric.color}
                        strokeWidth={3}
                        dot={{ fill: metric.color, r: 3 }}
                        activeDot={{ r: 5 }}
                        />
                    </LineChart>
                    </ResponsiveContainer>
                    <div className="mt-2 text-center">
                    <div className="text-2xl font-bold text-white">
                        {metricsHistory[metricsHistory.length - 1][metricKey].toFixed(2)}
                    </div>
                    <div className="text-xs text-purple-300">Current Value</div>
                    </div>
                </div>
                );
            })}
            </div>
        );
    }, [metricsHistory, visibleMetrics]);

    // Update metrics history whenever tree changes
    useEffect(() => {
        if (treeNodes.length > 0) {
            let metrics = calculateMetrics(treeNodes);
            setMetricsHistory(prev => {
                metrics.step = prev.length;
                return [...prev, { ...metrics }];
            });
        }
    }, [treeNodes.length]);

    return (
        <div className="mt-6 bg-white/10 backdrop-blur-lg rounded-2xl shadow-2xl p-6 border border-white/20">
            <div className="flex items-center justify-between mb-4">
                <h3 className="text-xl font-semibold text-white">Optimization Metrics</h3>
                
                {/* Metric Toggles */}
                <div className="flex gap-2">
                {Object.keys(metricDefinitions).map(key => (
                    <button
                    key={key}
                    onClick={() => setVisibleMetrics(prev => ({ ...prev, [key]: !prev[key as keyof VisibleMetrics] }))}
                    className={`px-3 py-1 rounded-lg text-sm font-medium transition-all ${
                        visibleMetrics[key as keyof VisibleMetrics]
                        ? 'bg-purple-500 text-white'
                        : 'bg-white/10 text-purple-300 hover:bg-white/20'
                    }`}
                    >
                    {metricDefinitions[key].label}
                    </button>
                ))}
                </div>
            </div>

            {metricsHistory.length > 0 ? (
                <>
                {/* Plots Grid - Memoized to prevent re-render on mouse move */}
                {metricsCharts}

                {/* Node Count Line */}
                {Object.values(visibleMetrics).filter(Boolean).length > 0 && (
                    <div className="mt-4 text-sm text-purple-300 text-center">
                    Total Molecules: {treeNodes.length}
                    </div>
                )}
                </>
            ) : (
                <div className="text-purple-300 text-center py-8">
                Waiting for metrics data...
                </div>
            )}
        </div>
    );
};
