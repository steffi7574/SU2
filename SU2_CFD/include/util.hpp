#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <stdarg.h>
#include <stdio.h>

/**
 * Conversion function which uses a format specifier for the string conversion.
 *
 * @param format        The format specifier like printf
 * @param list          The variable argument list for the format string
 *
 * @return  The output with the formated values
 */
std::string vformat(const char* format, va_list list);




/**
 * Conversion function which uses a format specifier for the string conversion.
 *
 * @param format        The format specifier like printf
 * @param ...           The values for the format string
 *
 * @return  The output with the formated values
 */
std::string format(const char* format, ...);
