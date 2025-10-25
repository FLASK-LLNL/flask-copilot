import { MarkdownTextProps } from "../types";

// Simple markdown-like parser for hover info
export const MarkdownText: React.FC<MarkdownTextProps> = ({ text }) => {
  const lines = text.split('\n');
  
  const parseInline = (line: string): string => {
    line = line.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');
    line = line.replace(/\*(.+?)\*/g, '<em>$1</em>');
    line = line.replace(/`(.+?)`/g, '<code class="bg-purple-900/50 px-1 rounded">$1</code>');
    return line;
  };
  
  return (
    <div className="space-y-2">
      {lines.map((line, idx) => {
        if (line.startsWith('# ')) {
          return <h3 key={idx} className="text-lg font-bold text-purple-200">{line.slice(2)}</h3>;
        } else if (line.startsWith('## ')) {
          return <h4 key={idx} className="text-base font-semibold text-purple-300">{line.slice(3)}</h4>;
        } else if (line.startsWith('- ')) {
          return (
            <li key={idx} className="ml-4 text-sm text-purple-100" dangerouslySetInnerHTML={{ __html: parseInline(line.slice(2)) }} />
          );
        } else if (line.trim() === '') {
          return <div key={idx} className="h-1" />;
        } else {
          return <p key={idx} className="text-sm text-purple-100" dangerouslySetInnerHTML={{ __html: parseInline(line) }} />;
        }
      })}
    </div>
  );
};
