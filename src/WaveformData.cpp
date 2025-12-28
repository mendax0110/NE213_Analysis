#include "../include/WaveformData.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

using namespace ne213;

std::vector<std::vector<double> > WaveformLoader::load_from_file(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file:" + filename);
    }

    std::vector<std::vector<double>> waveforms;
    std::string line;
    int line_count = 0;

    while (std::getline(file, line))
    {
        if (line_count % 3 == 0 && !line.empty())
        {
            auto data = parse_line(line);
            if (!data.empty())
            {
                waveforms.push_back(std::move(data));
            }
        }

        ++line_count;
    }

    return waveforms;
}

std::vector<double> WaveformLoader::parse_line(const std::string& line)
{
    std::string cleaned = line;

    if (cleaned.substr(0, 4) == "CH1:")
    {
        cleaned = cleaned.substr(4);
    }

    std::vector<double> data;
    std::stringstream ss(cleaned);
    std::string value;

    while (std::getline(ss, value, ','))
    {
        if (!value.empty())
        {
            try
            {
                double val = std::stod(value);
                val = (val > 128.0) ? (val - 255.0) : val;
                data.push_back(val);
            }
            catch (...)
            {
                continue;
            }
        }
    }

    if (data.size() > 350)
    {
        data.resize(350);
    }

    return data;
}
