#include "../include/WaveformData.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

using namespace ne213;

std::vector<std::vector<double>> WaveformLoader::load_from_file(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + filename);
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

FileType WaveformLoader::classify_file_type(const std::string& filename)
{
    std::string basename = filename;

    const size_t last_slash = filename.find_last_of("/\\");
    if (last_slash != std::string::npos)
    {
        basename = filename.substr(last_slash + 1);
    }

    std::string lower_basename;
    lower_basename.reserve(basename.size());
    for (const char c : basename)
    {
        lower_basename.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }

    if (lower_basename.find("mit") != std::string::npos)
    {
        return FileType::MIT;
    }
    if (lower_basename.find("ohne") != std::string::npos)
    {
        return FileType::OHNE;
    }

    return FileType::UNKNOWN;
}

std::string WaveformLoader::file_type_to_string(const FileType type)
{
    switch (type)
    {
        case FileType::MIT:
            return "MIT (with fusor)";
        case FileType::OHNE:
            return "OHNE (background)";
        case FileType::UNKNOWN:
        default:
            return "UNKNOWN";
    }
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

    if (data.size() > MAX_WAVEFORM_SAMPLES)
    {
        data.resize(MAX_WAVEFORM_SAMPLES);
    }

    return data;
}
