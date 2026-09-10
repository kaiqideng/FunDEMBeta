#include "datWriter.h"

#include <cctype>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace fundem
{
namespace
{

namespace fs = std::filesystem;

void createParentDirectory(const std::string& fileName)
{
    if (fileName.empty())
    {
        throw std::invalid_argument("DAT file name must not be empty.");
    }

    const fs::path parent = fs::path(fileName).parent_path();
    if (parent.empty())
    {
        return;
    }

    std::error_code error;
    fs::create_directories(parent, error);
    if (error)
    {
        throw std::runtime_error("Cannot create DAT output directory '" + parent.string() + "': " + error.message());
    }
}

bool containsWhitespace(const std::string& text)
{
    for (unsigned char character : text)
    {
        if (std::isspace(character) != 0)
        {
            return true;
        }
    }
    return false;
}

} // namespace

datWriter::datWriter(std::string fileName, std::vector<std::string> columnNames, bool appendExisting) : fileName_(std::move(fileName)), columnNames_(std::move(columnNames))
{
    if (fileName_.empty() || columnNames_.empty())
    {
        throw std::invalid_argument("DAT output requires a file name and at least one column.");
    }
    std::unordered_set<std::string> uniqueColumnNames;
    uniqueColumnNames.reserve(columnNames_.size());
    for (const std::string& columnName : columnNames_)
    {
        if (columnName.empty() || containsWhitespace(columnName) || !uniqueColumnNames.insert(columnName).second)
        {
            throw std::invalid_argument("DAT column names must be unique, non-empty, and contain no whitespace.");
        }
    }

    createParentDirectory(fileName_);
    bool writeNewHeader = true;
    if (appendExisting)
    {
        std::error_code error;
        writeNewHeader = !fs::exists(fileName_, error) || fs::file_size(fileName_, error) == 0;
        if (error)
        {
            throw std::runtime_error("Cannot inspect DAT file '" + fileName_ + "': " + error.message());
        }
        if (!writeNewHeader)
        {
            std::ifstream existing(fileName_);
            std::string header;
            std::getline(existing, header);
            std::ostringstream expectedHeader;
            expectedHeader << '#';
            for (const std::string& columnName : columnNames_)
            {
                expectedHeader << ' ' << columnName;
            }
            if (!existing && !existing.eof())
            {
                throw std::runtime_error("Cannot read DAT header: " + fileName_);
            }
            if (header != expectedHeader.str())
            {
                throw std::invalid_argument("Existing DAT header does not match the requested columns: " + fileName_);
            }
        }
    }

    const std::ios::openmode mode = std::ios::out | (appendExisting ? std::ios::app : std::ios::trunc);
    output_.open(fileName_, mode);
    if (!output_)
    {
        throw std::runtime_error("Cannot open DAT file for writing: " + fileName_);
    }
    output_ << std::scientific << std::setprecision(17);
    if (writeNewHeader)
    {
        writeHeader();
    }
}

void datWriter::writeHeader()
{
    output_ << '#';
    for (const std::string& columnName : columnNames_)
    {
        output_ << ' ' << columnName;
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing DAT header: " + fileName_);
    }
}

void datWriter::appendRow(const std::vector<Real>& values)
{
    if (values.size() != columnNames_.size())
    {
        throw std::invalid_argument("DAT row size does not match the column count.");
    }
    for (int columnIndex = 0; columnIndex < static_cast<int>(values.size()); ++columnIndex)
    {
        if (columnIndex > 0)
        {
            output_ << ' ';
        }
        output_ << values[columnIndex];
    }
    output_ << '\n';
    if (!output_)
    {
        throw std::runtime_error("Failed while writing DAT row: " + fileName_);
    }
}

void datWriter::appendRow(std::initializer_list<Real> values) { appendRow(std::vector<Real>(values)); }

void datWriter::flush()
{
    output_.flush();
    if (!output_)
    {
        throw std::runtime_error("Failed while flushing DAT file: " + fileName_);
    }
}

} // namespace fundem
