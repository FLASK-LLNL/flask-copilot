import { MarkdownTextProps } from "../types";

// Simple markdown-like parser for hover info
export const MarkdownText: React.FC<MarkdownTextProps> = ({ text }) => {
  const lines = text.split('\n');
  
  const parseInline = (line: string): string => {
    // Bold text
    line = line.replace(/\*\*(.+?)\*\*/g, '<strong class="markdown-strong">$1</strong>');
    // Italic text
    line = line.replace(/\*(.+?)\*/g, '<em class="markdown-em">$1</em>');
    // Inline code
    line = line.replace(/`(.+?)`/g, '<code class="markdown-code">$1</code>');
    return line;
  };
  
  return (
    <div className="markdown-content space-y-2">
      {lines.map((line, idx) => {
        // Heading level 1
        if (line.startsWith('# ')) {
          return (
            <h3 key={idx} className="markdown-heading-1">
              {line.slice(2)}
            </h3>
          );
        } 
        // Heading level 2
        else if (line.startsWith('## ')) {
          return (
            <h4 key={idx} className="markdown-heading-2">
              {line.slice(3)}
            </h4>
          );
        } 
        // List item
        else if (line.startsWith('- ')) {
          return (
            <li 
              key={idx} 
              className="markdown-list-item" 
              dangerouslySetInnerHTML={{ __html: parseInline(line.slice(2)) }} 
            />
          );
        } 
        // Empty line (spacer)
        else if (line.trim() === '') {
          return <div key={idx} className="markdown-spacer" />;
        } 
        // Regular paragraph
        else {
          return (
            <p 
              key={idx} 
              className="markdown-paragraph" 
              dangerouslySetInnerHTML={{ __html: parseInline(line) }} 
            />
          );
        }
      })}
    </div>
  );
};
