#pragma once

#include <Python.h> // has to be the first included header
#include <string>
#include <vector>

namespace StringUtility {

//! extract from a line in a file the literal that follows after a string key
//! and convert it to int
int getNumberAfterString(std::string line, std::string key);

//! remove beginning of line until (including) key, return removed substring
std::string extractUntil(std::string &line, std::string key);

//! remove whitespace (' ', '\t', '\n') at the beginning and end of the string
void trim(std::string &str);

//! output the values separated by spaces, after nValuesPerRow there will be a
//! line break, disabled if -1
template <typename IterType>
void outputValuesBlock(std::ostream &stream, IterType valuesBegin,
                       IterType valuesEnd, int nValuesPerRow = -1);

//! output the values separated by spaces, after nValuesPerRow there will be a
//! line break, disabled if -1, add the value 1 to each value
template <typename IterType>
void outputValuesBlockAdd1(std::ostream &stream, IterType valuesBegin,
                           IterType valuesEnd, int nValuesPerRow = -1);

//! replace the first occurence of from by to
std::string replace(std::string str, const std::string &from,
                    const std::string &to);

//! replace all occurences of from by to
std::string replaceAll(std::string str, const std::string &from,
                       const std::string &to);

//! for N=1 output <str>, for N=2 output <str>*<str>, for N=3 output
//! <str>*<str>*<str>
template <int N> std::string multiply(std::string str);

//! extract the basename of a file, i.e. remove leading path and trailing .*
std::string extractBasename(std::string str);

//! converts time to string object
std::string timeToString(const tm *const time);

//! return the human readable version of the result of typeid(<class>).name()
std::string demangle(const char *typeidName);

//! get the correct length of the string, also if it contains utf-8 based
//! unicode charactors, such as λ,γ etc.
std::size_t stringLength(std::string string);

} // namespace StringUtility

#include "utility/string_utility.tpp"
