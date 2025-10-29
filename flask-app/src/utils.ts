// General utility functions

export const copyToClipboard = async (text: string, fieldName: string, setCopiedField: (name: string | null) => void): Promise<void> => {
    try {
        await navigator.clipboard.writeText(text);
        setCopiedField(fieldName);
        setTimeout(() => setCopiedField(null), 2000);
    } catch (err) {
        console.error('Failed to copy:', err);
    }
};
