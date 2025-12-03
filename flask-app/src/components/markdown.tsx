import { MarkdownTextProps } from "../types";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { vscDarkPlus } from 'react-syntax-highlighter/dist/esm/styles/prism';

// Markdown renderer with custom styling and syntax highlighting
export const MarkdownText: React.FC<MarkdownTextProps> = ({ text }) => {
  return (
    <div className="markdown-content space-y-2">
      <ReactMarkdown
        remarkPlugins={[remarkGfm]}
        components={{
          h1: ({ children }) => <h3 className="markdown-heading-1">{children}</h3>,
          h2: ({ children }) => <h4 className="markdown-heading-2">{children}</h4>,
          h3: ({ children }) => <h5 className="markdown-heading-3">{children}</h5>,
          p: ({ children }) => <p className="markdown-paragraph">{children}</p>,
          ul: ({ children }) => <ul className="list-disc list-inside">{children}</ul>,
          ol: ({ children }) => <ol className="list-decimal list-inside">{children}</ol>,
          li: ({ children }) => <li className="markdown-list-item">{children}</li>,
          strong: ({ children }) => <strong className="markdown-strong">{children}</strong>,
          em: ({ children }) => <em className="markdown-em">{children}</em>,
          code: ({ node, inline, className, children, ...props }: any) => {
            const match = /language-(\w+)/.exec(className || '');

            // Multi-line code block with syntax highlighting
            if (!inline && match) {
              return (
                <SyntaxHighlighter
                  style={vscDarkPlus}
                  language={match[1]}
                  PreTag="div"
                  {...props}
                >
                  {String(children).replace(/\n$/, '')}
                </SyntaxHighlighter>
              );
            }

            // Single- or Multi-line code block without language
            return (
              <code {...props}>{children}</code>
            );

          },
          table: ({ children }) => <table className="markdown-table">{children}</table>,
          th: ({ children }) => <th className="markdown-table-header">{children}</th>,
          td: ({ children }) => <td className="markdown-table-cell">{children}</td>,
        }}
      >
        {text}
      </ReactMarkdown>
    </div>
  );
};
