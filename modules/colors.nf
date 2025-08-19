/*
 * Color utilities for NextITS pipeline
 * Provides ANSI color codes and helper functions for terminal output
 */

/**
 * Get color codes map based on monochrome setting
 * @param monochrome Whether to disable colors (from params.monochrome_logs)
 * @return Map of color names to ANSI codes
 */
def getColors(boolean monochrome = false) {
    return monochrome ? [:] : [
        // Basic colors
        black:   "\033[0;30m",
        red:     "\033[0;31m",
        green:   "\033[0;32m",
        yellow:  "\033[0;33m",
        blue:    "\033[0;34m",
        purple:  "\033[0;35m",
        cyan:    "\033[0;36m",
        white:   "\033[0;37m",
        
        // Bright colors
        bright_black:   "\033[0;90m",
        bright_red:     "\033[0;91m",
        bright_green:   "\033[0;92m",
        bright_yellow:  "\033[0;93m",
        bright_blue:    "\033[0;94m",
        bright_purple:  "\033[0;95m",
        bright_cyan:    "\033[0;96m",
        bright_white:   "\033[0;97m",
        
        // Text formatting
        bold:      "\033[1m",
        dim:       "\033[2m",
        italic:    "\033[3m",
        underline: "\033[4m",
        blink:     "\033[5m",
        reverse:   "\033[7m",
        
        // Reset
        reset:     "\033[0m"
    ]
}

/**
 * Apply color formatting to text
 * @param text The text to colorize
 * @param color The color name (e.g., 'red', 'green', 'bold')
 * @param monochrome Whether to disable colors
 * @return Formatted text string
 */
def colorize(String text, String color, boolean monochrome = false) {
    def colors = getColors(monochrome)
    if (!colors[color]) {
        return text
    }
    return "${colors[color]}${text}${colors.reset}"
}

/**
 * Apply multiple color/format combinations to text
 * @param text The text to colorize
 * @param formats List of format names (e.g., ['red', 'bold'])
 * @param monochrome Whether to disable colors
 * @return Formatted text string
 */
def colorizeMultiple(String text, List<String> formats, boolean monochrome = false) {
    def colors = getColors(monochrome)
    if (monochrome || !formats) {
        return text
    }
    
    def prefix = formats.findAll { colors[it] }.collect { colors[it] }.join('')
    return "${prefix}${text}${colors.reset}"
}

/**
 * Create an error message with red coloring
 * @param message The error message text
 * @param monochrome Whether to disable colors
 * @return Formatted error message
 */
def errorMsg(String message, boolean monochrome = false) {
    return colorizeMultiple("ERROR: ${message}", ['red', 'bold'], monochrome)
}

/**
 * Create a warning message with yellow coloring
 * @param message The warning message text
 * @param monochrome Whether to disable colors
 * @return Formatted warning message
 */
def warningMsg(String message, boolean monochrome = false) {
    return colorizeMultiple("WARNING: ${message}", ['yellow', 'bold'], monochrome)
}

/**
 * Create an info message with cyan coloring
 * @param message The info message text
 * @param monochrome Whether to disable colors
 * @return Formatted info message
 */
def infoMsg(String message, boolean monochrome = false) {
    return colorize(message, 'cyan', monochrome)
}

/**
 * Create a success message with green coloring
 * @param message The success message text
 * @param monochrome Whether to disable colors
 * @return Formatted success message
 */
def successMsg(String message, boolean monochrome = false) {
    return colorizeMultiple(message, ['green', 'bold'], monochrome)
}
