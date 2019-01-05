#pragma once
#ifndef pathutil_hpp
#define pathutil_hpp

void mkpath(const char *prefix, const char *filename, char *path);

void mkrankpath(int rank, const char *prefix, const char *filename, char *path);

#endif // pathutil_hpp
