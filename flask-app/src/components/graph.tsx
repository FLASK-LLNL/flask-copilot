// Molecule graph view
import { Loader2, Move } from "lucide-react";
import { useRef, useCallback, useEffect, useState } from "react";
import { BOX_HEIGHT, BOX_WIDTH, NODE_STYLES } from "../constants";
import {MoleculeGraphProps, MoleculeGraphState, Position, TreeNode} from "../types";
import { MoleculeSVG } from "./molecule";
import { MarkdownText } from "./markdown";

export const useGraphState = (): MoleculeGraphState => {
    const [offset, setOffset] = useState<Position>({ x: 50, y: 50 });
    const [zoom, setZoom] = useState<number>(1);
  
    return { offset, setOffset, zoom, setZoom };
};

export const MoleculeGraph: React.FC<MoleculeGraphProps> = ({nodes, edges, ctx, autoZoom, setAutoZoom, handleNodeClick, rdkitModule, offset, setOffset, zoom, setZoom}) => {
    const [isDragging, setIsDragging] = useState<boolean>(false);
    const [dragStart, setDragStart] = useState<Position>({ x: 0, y: 0 });
    const [hoveredNode, setHoveredNode] = useState<TreeNode | null>(null);
    const [mousePos, setMousePos] = useState<Position>({ x: 0, y: 0 });
    const containerRef = useRef<HTMLDivElement>(null);

    const getNode = useCallback((nodeId: string): TreeNode | undefined => {
        return nodes.find(n => n.id === nodeId);
    }, [nodes]);


    const handleMouseDown = (e: React.MouseEvent<HTMLDivElement>): void => {
    if (e.button !== 0) return;
    setIsDragging(true);
    setDragStart({ x: e.clientX - offset.x, y: e.clientY - offset.y });
    e.preventDefault();
    };

    const handleMouseMove = (e: MouseEvent): void => {
    if (!isDragging) return;
    setOffset({
        x: e.clientX - dragStart.x,
        y: e.clientY - dragStart.y
    });
    };

    const handleMouseUp = (): void => {
    setIsDragging(false);
    };

    const handleWheel = (e: WheelEvent): void => {
    e.preventDefault();
    
    if (autoZoom) {
        setAutoZoom(false);
    }
    
    const rect = containerRef.current?.getBoundingClientRect();
        if (!rect) return;
        
        const mouseX = e.clientX - rect.left;
        const mouseY = e.clientY - rect.top;
        
        const delta = -e.deltaY * 0.001;
        const newZoom = Math.min(Math.max(0.1, zoom + delta), 3);
        
        const zoomRatio = newZoom / zoom;
        const newOffsetX = mouseX - (mouseX - offset.x) * zoomRatio;
        const newOffsetY = mouseY - (mouseY - offset.y) * zoomRatio;
        
        setZoom(newZoom);
        setOffset({ x: newOffsetX, y: newOffsetY });
    };

    const resetZoom = (): void => {
        setZoom(1);
        setOffset({ x: 50, y: 50 });
    };

    const fitToView = (): void => {
    if (!containerRef.current || nodes.length === 0) return;
    
    const padding = 80;
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    
    const nodeElements = containerRef.current.querySelectorAll('[data-node-id]');
    
    if (nodeElements.length === 0) {
        nodes.forEach(node => {
        minX = Math.min(minX, node.x);
        maxX = Math.max(maxX, node.x + 140);
        minY = Math.min(minY, node.y);
        maxY = Math.max(maxY, node.y + 140);
        });
    } else {
        nodeElements.forEach((element) => {
        const nodeId = element.getAttribute('data-node-id');
        const node = nodes.find(n => n.id === nodeId);
        
        if (node) {
            const actualWidth = (element as HTMLElement).offsetWidth;
            const actualHeight = (element as HTMLElement).offsetHeight;
            
            minX = Math.min(minX, node.x);
            maxX = Math.max(maxX, node.x + actualWidth);
            minY = Math.min(minY, node.y);
            maxY = Math.max(maxY, node.y + actualHeight);
        }
        });
    }
    
    const contentWidth = maxX - minX + padding * 2;
    const contentHeight = maxY - minY + padding * 2;
    
    const rect = containerRef.current.getBoundingClientRect();
    const viewWidth = rect.width;
    const viewHeight = rect.height;
    
    const scaleX = viewWidth / contentWidth;
    const scaleY = viewHeight / contentHeight;
    const newZoom = Math.min(scaleX, scaleY, 1);
    
    const contentCenterX = (minX + maxX) / 2;
    const contentCenterY = (minY + maxY) / 2;
    
    const newOffsetX = viewWidth / 2 - contentCenterX * newZoom;
    const newOffsetY = viewHeight / 2 - contentCenterY * newZoom;
    
    setZoom(newZoom);
    setOffset({ x: newOffsetX, y: newOffsetY });
    };

    useEffect(() => {
    if (isDragging) {
        window.addEventListener('mousemove', handleMouseMove);
        window.addEventListener('mouseup', handleMouseUp);
        return () => {
        window.removeEventListener('mousemove', handleMouseMove);
        window.removeEventListener('mouseup', handleMouseUp);
        };
    }
    }, [isDragging, dragStart]);

    // Add wheel listener with passive: false to allow preventDefault
    useEffect(() => {
    const container = containerRef.current;
    if (container) {
        container.addEventListener('wheel', handleWheel, { passive: false });
        return () => {
        container.removeEventListener('wheel', handleWheel);
        };
    }
    }, [autoZoom, zoom, offset]);

    useEffect(() => {
    if (autoZoom && nodes.length > 0) {
        requestAnimationFrame(() => {
            requestAnimationFrame(() => {
                fitToView();
            });
        });
    }
    }, [nodes, autoZoom]);


    const getCurvedPath = (from: TreeNode | undefined, to: TreeNode | undefined): string => {
        if (!from || !to) return '';

        const startX = from.x + BOX_WIDTH;
        const startY = from.y + BOX_HEIGHT / 2;
        const endX = to.x;
        const endY = to.y + BOX_HEIGHT / 2;

        const controlX1 = startX + (endX - startX) * 0.3;
        const controlX2 = startX + (endX - startX) * 0.7;

        return `M ${startX} ${startY} C ${controlX1} ${startY}, ${controlX2} ${endY}, ${endX} ${endY}`;
    };

    const getCurveMidpoint = (from: TreeNode | undefined, to: TreeNode | undefined): Position => {
        if (!from || !to) return { x: 0, y: 0 };
        const startX = from.x + BOX_WIDTH;
        const startY = from.y + BOX_HEIGHT / 2;
        const endX = to.x;
        const endY = to.y + BOX_HEIGHT / 2;

        const t = 0.5;
        const controlX1 = startX + (endX - startX) * 0.3;
        const controlX2 = startX + (endX - startX) * 0.7;

        const x = Math.pow(1-t, 3) * startX + 
                    3 * Math.pow(1-t, 2) * t * controlX1 + 
                    3 * (1-t) * Math.pow(t, 2) * controlX2 + 
                    Math.pow(t, 3) * endX;

        const y = Math.pow(1-t, 3) * startY + 
                    3 * Math.pow(1-t, 2) * t * startY + 
                    3 * (1-t) * Math.pow(t, 2) * endY + 
                    Math.pow(t, 3) * endY;

        return { x, y };
    };

    return (
        <>
            <div 
                ref={containerRef}
                className={`relative w-full h-full overflow-hidden ${
                autoZoom ? 'cursor-default' : isDragging ? 'cursor-grabbing' : 'cursor-grab'
                }`}
                onMouseDown={handleMouseDown}
                style={{ userSelect: 'none' }}
            >
                <div
                className="absolute"
                style={{
                    transform: `translate(${offset.x}px, ${offset.y}px) scale(${zoom})`,
                    transformOrigin: '0 0',
                    transition: isDragging ? 'none' : 'transform 0.1s ease-out',
                    width: '3000px',
                    height: '2000px'
                }}
                >
                {edges.filter(edge => getNode(edge.fromNode) && getNode(edge.toNode)).map((edge, idx) => {
                    const midpoint = getCurveMidpoint(getNode(edge.fromNode), getNode(edge.toNode));
                    return (
                    <div key={edge.id} className="absolute pointer-events-none">
                        <svg className="absolute" style={{ width: '3000px', height: '2000px', top: 0, left: 0 }}>
                        <g className="animate-fadeIn" style={{ animationDelay: `${idx * 50}ms` }}>
                            <path
                            d={getCurvedPath(getNode(edge.fromNode), getNode(edge.toNode))}
                            stroke={edge.status === 'computing' ? '#F59E0B' : '#8B5CF6'}
                            strokeWidth="3"
                            fill="none"
                            strokeDasharray={edge.status === 'computing' ? '5,5' : 'none'}
                            className={edge.status === 'computing' ? 'animate-dash' : ''}
                            opacity="0.8"
                            />
                            <circle cx={getNode(edge.toNode)!.x + 10} cy={getNode(edge.toNode)!.y + 50} r="5" fill={edge.status === 'computing' ? '#F59E0B' : '#EC4899'} />
                        </g>
                        </svg>
                        
                        <div className="absolute pointer-events-auto" style={{ left: `${midpoint.x}px`, top: `${midpoint.y}px`, transform: 'translate(-50%, -50%)' }}>
                        {edge.label && (
                            <div className={`px-3 py-1.5 rounded-lg text-xs font-semibold shadow-lg ${edge.status === 'computing' ? 'bg-amber-500 text-white animate-pulse' : 'bg-purple-500 text-white'}`}>
                            {edge.status === 'computing' && <Loader2 className="w-3 h-3 inline mr-1 animate-spin" />}
                            {edge.label}
                            </div>
                        )}
                        </div>
                    </div>
                    );
                })}

                {nodes.map((node, idx) => (
                    <div
                    key={node.id}
                    data-node-id={node.id}
                    className="absolute animate-fadeInScale"
                    style={{ 
                        left: `${node.x}px`, 
                        top: `${node.y}px`, 
                        width: `${BOX_WIDTH*1.05}px`,
                        animationDelay: `${idx * 100}ms`,
                        transition: 'left 0.6s cubic-bezier(0.34, 1.56, 0.64, 1), top 0.6s cubic-bezier(0.34, 1.56, 0.64, 1)'
                    }}
                    onMouseEnter={(e) => {
                        setHoveredNode(node);
                        setMousePos({ x: e.clientX, y: e.clientY });
                    }}
                    onMouseMove={(e) => {
                        if (hoveredNode?.id === node.id) {
                        setMousePos({ x: e.clientX, y: e.clientY });
                        }
                    }}
                    onMouseLeave={() => setHoveredNode(null)}
                    onClick={(e) => handleNodeClick(e, node)}
                    >
                    <div className={`bg-gradient-to-br backdrop-blur-sm rounded-xl p-3 border-2 shadow-lg hover:shadow-2xl hover:scale-105 transition-all pointer-events-auto cursor-pointer ${NODE_STYLES[node.highlight || 'normal']}`} style={{ textAlign: 'center' }}>
                        <MoleculeSVG smiles={node.smiles} height={80} rdkitModule={rdkitModule} />
                        <div className="mt-2 text-center">
                        <div className="text-xs font-semibold text-purple-200 bg-black/30 rounded px-2 py-1 whitespace-pre-line">{node.label}</div>
                        </div>
                    </div>
                    </div>
                ))}
                </div>

                {nodes.length > 0 && !isDragging && (
                <div className="absolute bottom-4 left-4 space-y-2 pointer-events-none">
                    <div className="bg-black/50 backdrop-blur-sm rounded-lg px-3 py-2 text-purple-200 text-sm flex items-center gap-2">
                    <Move className="w-4 h-4" />
                    {autoZoom ? 'Auto-zoom active • Drag to pan' : 'Drag to pan • Scroll to zoom'}
                    </div>
                    <div className="bg-black/50 backdrop-blur-sm rounded-lg px-3 py-2 text-purple-200 text-sm flex items-center justify-between gap-3 pointer-events-auto">
                    <span>Zoom: {(zoom * 100).toFixed(0)}%</span>
                    {!autoZoom && (
                        <button
                        onClick={(e) => { e.stopPropagation(); resetZoom(); }}
                        onMouseDown={(e) => e.stopPropagation()}
                        className="text-xs bg-purple-500 hover:bg-purple-600 px-2 py-1 rounded transition-colors"
                        >
                        Reset
                        </button>
                    )}
                    </div>
                </div>
                )}
            </div>
            {hoveredNode && !ctx && (
                <div className="fixed z-50 pointer-events-none" style={{ left: `${mousePos.x + 20}px`, top: `${mousePos.y + 20}px`, maxWidth: '400px' }}>
                <div className="bg-gradient-to-br from-slate-800 to-purple-900 border-2 border-purple-400 rounded-xl shadow-2xl p-4 max-h-96 overflow-y-auto">
                    { /* <div className="text-xs text-purple-400 mb-2">Debug: x={mousePos.x}, y={mousePos.y}</div> */ }
                    <MarkdownText text={hoveredNode.hoverInfo} />
                </div>
                </div>
            )}
        </>
    );
};
