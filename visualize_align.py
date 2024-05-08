
def wrap_text(text, max_line_length=60):
    """Wraps a long text into lines with a maximum character limit.

    Args:
        text: The text to be wrapped.
        max_line_length: The maximum number of characters allowed per line.

    Returns:
        A list of strings, where each string is a line of the wrapped text.
    """

    lines = []
    current_line = ""

    for char in text:
        # If adding the character would exceed the line length, add the current line to the list and start a new one
        if len(current_line) + len(char) > max_line_length:
            lines.append(current_line)
            current_line = ""
  
        # Add the character to the current line
        current_line += char

    # Add the last line (if it's not empty)
    if current_line:
        lines.append(current_line)

    return lines

if __name__ == "__main__":
    pass