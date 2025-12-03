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
    sascore: { 
        label: 'Synthesizability Score', 
        color: '#48adff',
        calculate: (nodes: TreeNode[]) => {
            // Get the lowest SA score among all nodes
            const scores = nodes.map(n => n.sascore || 0);
            return scores.length > 0 ? Math.min(...scores) : 0;
        }
    },
    bandgap: { 
        label: 'Band Gap (eV)', 
        color: '#F59E0B',
        calculate: (nodes: TreeNode[]) => {
            const bandgaps = nodes.map(n => n.bandgap || 0).filter(b => b > 0);
            return bandgaps.length > 0 ? Math.max(...bandgaps) : 0;
        }
    },
    density: {
        label: 'Molecular Density (g/cm³)',
        color: '#10B981',
        calculate: (nodes: TreeNode[]) => {
            // Get the highest density among all nodes
            const densities = nodes.map(n => n.density || 0).filter(d => d > 0);
            return densities.length > 0 ? Math.max(...densities) : 0;
        }
    },
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
        sascore: false,
        density: false,
        yield: false,
    });

    return {metricsHistory, setMetricsHistory, visibleMetrics, setVisibleMetrics};
};

export const MetricsDashboard: React.FC<MetricsDashboardProps> = ({treeNodes, metricsHistory, setMetricsHistory, visibleMetrics, setVisibleMetrics}) => {
    // Memoize the metrics charts to prevent re-render on mouse move
    const metricsCharts = useMemo(() => {
        if (metricsHistory.length === 0) return null;

        const visibleCount = Object.values(visibleMetrics).filter(Boolean).length;
        
        return (
            <div className={`metrics-grid ${visibleCount > 1 ? 'metrics-grid-multi' : 'metrics-grid-single'}`}>
            {Object.keys(metricDefinitions).filter(key => visibleMetrics[key as keyof VisibleMetrics]).map(metricKey => {
                const metric = metricDefinitions[metricKey];
                return (
                <div key={metricKey} className="metric-card">
                    <h4 className="metric-title">{metric.label}</h4>
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
                    <div className="metric-value">
                    <div className="metric-value-number">
                        {metricsHistory[metricsHistory.length - 1][metricKey].toFixed(2)}
                    </div>
                    <div className="metric-value-label">Current Value</div>
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
            const metrics = calculateMetrics(treeNodes);
            setMetricsHistory(prev => {
                // Only add if this represents a new state
                const lastMetrics = prev[prev.length - 1];

                // Check if metrics actually changed
                const metricsChanged = !lastMetrics ||
                    Object.keys(metricDefinitions).some(key =>
                        Math.abs((lastMetrics[key] || 0) - (metrics[key] || 0)) > 0.0001
                    ) ||
                    lastMetrics.nodeCount !== metrics.nodeCount;

                if (!metricsChanged) {
                    return prev;
                }
                metrics.step = prev.length;
                return [...prev, { ...metrics }];
            });
        }
    }, [treeNodes]);

    return (
        <div className="card card-padding mt-6">
            <div className="flex-between mb-4">
                <h3 className="heading-3">Optimization Metrics</h3>

                {/* Metric Toggles */}
                <div className="flex gap-2">
                {Object.keys(metricDefinitions).map(key => (
                    <button
                    key={key}
                    onClick={() => setVisibleMetrics(prev => ({ ...prev, [key]: !prev[key as keyof VisibleMetrics] }))}
                    className={`btn btn-sm ${
                        visibleMetrics[key as keyof VisibleMetrics]
                        ? 'bg-primary text-primary'
                        : 'btn-tertiary text-tertiary'
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
                    <div className="mt-4 text-sm text-tertiary text-center">
                    Total Molecules: {treeNodes.length}
                    </div>
                )}
                </>
            ) : (
                <div className="text-tertiary text-center py-8">
                Waiting for metrics data...
                </div>
            )}
        </div>
    );
};
