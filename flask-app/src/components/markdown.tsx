import { MarkdownTextProps } from "../types";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { vscDarkPlus } from 'react-syntax-highlighter/dist/esm/styles/prism';

// Define components once outside the component to prevent recreation
const markdownComponents = {
  h1: ({ children }: any) => <h3 className="markdown-heading-1">{children}</h3>,
  h2: ({ children }: any) => <h4 className="markdown-heading-2">{children}</h4>,
  h3: ({ children }: any) => <h5 className="markdown-heading-3">{children}</h5>,
  p: ({ children }: any) => <p className="markdown-paragraph">{children}</p>,
  ul: ({ children }: any) => <ul className="list-disc list-inside">{children}</ul>,
  ol: ({ children }: any) => <ol className="list-decimal list-inside">{children}</ol>,
  li: ({ children }: any) => <li className="markdown-list-item">{children}</li>,
  strong: ({ children }: any) => <strong className="markdown-strong">{children}</strong>,
  em: ({ children }: any) => <em className="markdown-em">{children}</em>,
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
  table: ({ children }: any) => <table className="markdown-table">{children}</table>,
  th: ({ children }: any) => <th className="markdown-table-header">{children}</th>,
  td: ({ children }: any) => <td className="markdown-table-cell">{children}</td>,
};

// Markdown renderer with custom styling and syntax highlighting
export const MarkdownText: React.FC<MarkdownTextProps> = ({ text }) => {
  return (
    <div className="markdown-content space-y-2">
      <ReactMarkdown
        remarkPlugins={[remarkGfm]}
        components={markdownComponents}
      >
        {text}
      </ReactMarkdown>
    </div>
  );
};
