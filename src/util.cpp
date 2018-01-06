//
// Created by adyck on 1/6/18.
//
#include <cstdlib>
#include <string>

void stripColons(std::string &x)
{
    auto it = std::remove_if(std::begin(x),std::end(x),[](char c){return (c == ':');});
    x.erase(it, std::end(x));
}
