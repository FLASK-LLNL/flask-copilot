import { MarkdownTextProps } from "../types";

// Simple markdown-like parser
export const MarkdownText: React.FC<MarkdownTextProps> = ({ text }) => {
  const parseInline = (line: string): string => {
    // Bold text
    line = line.replace(/\*\*(.+?)\*\*/g, '<strong class="markdown-strong">$1</strong>');
    // Italic text
    line = line.replace(/\*(.+?)\*/g, '<em class="markdown-em">$1</em>');
    // Inline code
    line = line.replace(/`(.+?)`/g, '<code class="markdown-code">$1</code>');
    return line;
  };

  const parseBlocks = () => {
    const lines = text.split('\n');
    const blocks: JSX.Element[] = [];
    let i = 0;

    while (i < lines.length) {
      const line = lines[i];

      // Multiline code block
      if (line.startsWith('```')) {
        const codeLines: string[] = [];
        i++; // Skip opening ```
        while (i < lines.length && !lines[i].startsWith('```')) {
          codeLines.push(lines[i]);
          i++;
        }
        blocks.push(
          <pre key={blocks.length} className="markdown-code-block">
            <code>{codeLines.join('\n')}</code>
          </pre>
        );
        i++; // Skip closing ```
        continue;
      }

      // Table detection (line contains pipes)
      if (line.includes('|') && line.trim().startsWith('|')) {
        const tableLines: string[] = [line];
        let j = i + 1;
        // Collect consecutive table lines
        while (j < lines.length && lines[j].includes('|') && lines[j].trim().startsWith('|')) {
          tableLines.push(lines[j]);
          j++;
        }
        
        if (tableLines.length >= 2) {
          // Parse table
          const rows = tableLines.map(l => 
            l.split('|').slice(1, -1).map(cell => cell.trim())
          );
          const headers = rows[0];
          const dataRows = rows.slice(2); // Skip separator row
          
          blocks.push(
            <table key={blocks.length} className="markdown-table">
              <thead>
                <tr>
                  {headers.map((header, idx) => (
                    <th key={idx} className="markdown-table-header">
                      {header}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {dataRows.map((row, rowIdx) => (
                  <tr key={rowIdx}>
                    {row.map((cell, cellIdx) => (
                      <td 
                        key={cellIdx} 
                        className="markdown-table-cell"
                        dangerouslySetInnerHTML={{ __html: parseInline(cell) }}
                      />
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          );
          i = j;
          continue;
        }
      }

      // Heading level 1
      if (line.startsWith('# ')) {
        blocks.push(
          <h3 key={blocks.length} className="markdown-heading-1">
            {line.slice(2)}
          </h3>
        );
      } 
      // Heading level 2
      else if (line.startsWith('## ')) {
        blocks.push(
          <h4 key={blocks.length} className="markdown-heading-2">
            {line.slice(3)}
          </h4>
        );
      }
      // Heading level 3
      else if (line.startsWith('### ')) {
        blocks.push(
          <h5 key={blocks.length} className="markdown-heading-3">
            {line.slice(4)}
          </h5>
        );
      }
      // Numbered list item
      else if (/^\d+\.\s/.test(line)) {
        const match = line.match(/^\d+\.\s(.+)/);
        blocks.push(
          <li 
            key={blocks.length} 
            className="markdown-numbered-list-item" 
            dangerouslySetInnerHTML={{ __html: parseInline(match![1]) }} 
          />
        );
      }
      // Bullet list item (- or *)
      else if (line.startsWith('- ') || line.startsWith('* ')) {
        blocks.push(
          <li 
            key={blocks.length} 
            className="markdown-list-item" 
            dangerouslySetInnerHTML={{ __html: parseInline(line.slice(2)) }} 
          />
        );
      } 
      // Empty line (spacer)
      else if (line.trim() === '') {
        blocks.push(<div key={blocks.length} className="markdown-spacer" />);
      } 
      // Regular paragraph
      else {
        blocks.push(
          <p 
            key={blocks.length} 
            className="markdown-paragraph" 
            dangerouslySetInnerHTML={{ __html: parseInline(line) }} 
          />
        );
      }
      
      i++;
    }

    return blocks;
  };
  
  return (
    <div className="markdown-content space-y-2">
      {parseBlocks()}
    </div>
  );
};
