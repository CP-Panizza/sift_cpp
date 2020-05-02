//
// Created by Administrator on 2020/4/22.
//

#ifndef EL_UTILS_H
#define EL_UTILS_H



#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>


std::vector<std::string> split(std::string str, std::string pattern);

bool contain(std::string str, std::string target);

bool file_exists(const std::string &name);

bool dir_exists(std::string path);

long file_size(const char *filepath);

void trim_space(std::string &s);

std::string &replace_all(std::string &str, const std::string &old_value, const std::string &new_value);


#endif //EL_UTILS_H
