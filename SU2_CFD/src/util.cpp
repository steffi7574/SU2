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
std::string vformat(const char* format, va_list list) {
    const int bufferSize = 200;
    char buffer[bufferSize];

    // copy the list if we need to iterate through the variables again
    va_list listCpy;
    va_copy(listCpy, list);


    int outSize = vsnprintf(buffer, bufferSize, format, list);

    std::string result;
    if(outSize + 1 > bufferSize) {
        char* newBuffer = new char[outSize + 1];

        outSize = vsnprintf(newBuffer, outSize + 1, format, listCpy);

        result = newBuffer;

        delete [] newBuffer;
    } else {
        result = buffer;
    }

    // cleanup the copied list
    va_end (listCpy);

    return result;
}

/**
 * Conversion function which uses a format specifier for the string conversion.
 *
 * @param format        The format specifier like printf
 * @param ...           The values for the format string
 *
 * @return  The output with the formated values
 */
std::string format(const char* format, ...) {
    va_list list;
    va_start(list, format);
    std::string output = vformat(format, list);
    va_end(list);

    return output;
}
