#include "vtuWriter.h"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string_view>
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
        throw std::invalid_argument("VTU file name must not be empty.");
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
        throw std::runtime_error("Cannot create VTU output directory '" + parent.string() + "': " + error.message());
    }
}

std::string xmlEscape(std::string_view text)
{
    std::string escaped;
    escaped.reserve(text.size());
    for (char character : text)
    {
        switch (character)
        {
        case '&':
            escaped += "&amp;";
            break;
        case '<':
            escaped += "&lt;";
            break;
        case '>':
            escaped += "&gt;";
            break;
        case '"':
            escaped += "&quot;";
            break;
        case '\'':
            escaped += "&apos;";
            break;
        default:
            escaped += character;
            break;
        }
    }
    return escaped;
}

template <class Value> void writeASCIIValues(std::ostream& output, const std::vector<std::uint8_t>& bytes)
{
    if (bytes.size() % sizeof(Value) != 0)
    {
        throw std::logic_error("VTU data byte count does not match its declared type.");
    }
    if constexpr (std::is_floating_point_v<Value>)
    {
        output << std::setprecision(std::numeric_limits<Value>::max_digits10);
    }
    const std::size_t valueCount = bytes.size() / sizeof(Value);
    for (std::size_t index = 0; index < valueCount; ++index)
    {
        Value value;
        std::memcpy(&value, bytes.data() + index * sizeof(Value), sizeof(Value));
        if (index > 0)
        {
            output.put(' ');
        }
        if constexpr (sizeof(Value) == 1 && std::is_integral_v<Value>)
        {
            output << static_cast<int>(value);
        }
        else
        {
            output << value;
        }
    }
}

std::vector<vtuWriter::Real> flattenVectors(const std::vector<vtuWriter::Vec3>& vectors)
{
    std::vector<vtuWriter::Real> flattened;
    flattened.reserve(3 * vectors.size());
    for (const vtuWriter::Vec3& vector : vectors)
    {
        flattened.push_back(vector.x);
        flattened.push_back(vector.y);
        flattened.push_back(vector.z);
    }
    return flattened;
}

} // namespace

vtuWriter::vtuWriter(std::string fileName) : fileName_(std::move(fileName))
{
    if (fileName_.empty())
    {
        throw std::invalid_argument("VTU file name must not be empty.");
    }
}

vtuWriter::dataArray vtuWriter::makeDataArray(const char* type, const std::string& name, int componentCount, const void* data, std::size_t byteCount)
{
    if (type == nullptr || componentCount <= 0 || (byteCount > 0 && data == nullptr))
    {
        throw std::invalid_argument("Invalid VTU data array.");
    }

    dataArray array;
    array.type_ = type;
    array.name_ = name;
    array.componentCount_ = componentCount;
    array.bytes_.resize(byteCount);
    if (byteCount > 0)
    {
        std::memcpy(array.bytes_.data(), data, byteCount);
    }
    return array;
}

void vtuWriter::writeDataArray(std::ostream& output, const dataArray& array, vtuFormat format, std::uint64_t& appendedOffset)
{
    output << "<DataArray type=\"" << array.type_ << "\"";
    if (!array.name_.empty())
    {
        output << " Name=\"" << xmlEscape(array.name_) << "\"";
    }
    if (array.componentCount_ > 1)
    {
        output << " NumberOfComponents=\"" << array.componentCount_ << "\"";
    }
    if (format == vtuFormat::binary)
    {
        output << " format=\"appended\" offset=\"" << appendedOffset << "\"/>\n";
        appendedOffset += sizeof(std::uint64_t) + static_cast<std::uint64_t>(array.bytes_.size());
    }
    else
    {
        output << " format=\"ascii\">\n";
        writeASCIIData(output, array);
        output << "\n</DataArray>\n";
    }
}

void vtuWriter::writeASCIIData(std::ostream& output, const dataArray& array)
{
    if (array.type_ == "Float32")
    {
        writeASCIIValues<float>(output, array.bytes_);
    }
    else if (array.type_ == "Float64")
    {
        writeASCIIValues<double>(output, array.bytes_);
    }
    else if (array.type_ == "Int8")
    {
        writeASCIIValues<std::int8_t>(output, array.bytes_);
    }
    else if (array.type_ == "UInt8")
    {
        writeASCIIValues<std::uint8_t>(output, array.bytes_);
    }
    else if (array.type_ == "Int16")
    {
        writeASCIIValues<std::int16_t>(output, array.bytes_);
    }
    else if (array.type_ == "UInt16")
    {
        writeASCIIValues<std::uint16_t>(output, array.bytes_);
    }
    else if (array.type_ == "Int32")
    {
        writeASCIIValues<std::int32_t>(output, array.bytes_);
    }
    else if (array.type_ == "UInt32")
    {
        writeASCIIValues<std::uint32_t>(output, array.bytes_);
    }
    else if (array.type_ == "Int64")
    {
        writeASCIIValues<std::int64_t>(output, array.bytes_);
    }
    else if (array.type_ == "UInt64")
    {
        writeASCIIValues<std::uint64_t>(output, array.bytes_);
    }
    else
    {
        throw std::logic_error("Unsupported VTU data type: " + array.type_);
    }
}

void vtuWriter::writeAppendedData(std::ostream& output, const dataArray& array)
{
    const std::uint64_t byteCount = static_cast<std::uint64_t>(array.bytes_.size());
    output.write(reinterpret_cast<const char*>(&byteCount), sizeof(byteCount));
    if (!array.bytes_.empty())
    {
        output.write(reinterpret_cast<const char*>(array.bytes_.data()), static_cast<std::streamsize>(array.bytes_.size()));
    }
}

void vtuWriter::appendDataArray(std::vector<dataArray>& arrays, const char* type, const std::string& name, int componentCount, const void* data, std::size_t byteCount, int valueCount, int entityCount)
{
    if (name.empty() || componentCount <= 0 || valueCount != entityCount * componentCount)
    {
        throw std::invalid_argument("VTU field size does not match its point or cell count.");
    }
    if (std::any_of(arrays.begin(), arrays.end(), [&](const dataArray& array) { return array.name_ == name; }))
    {
        throw std::invalid_argument("Duplicate VTU field name: " + name);
    }
    arrays.push_back(makeDataArray(type, name, componentCount, data, byteCount));
}

void vtuWriter::setPoints(const std::vector<Vec3>& points)
{
    const int newPointCount = static_cast<int>(points.size());
    if ((cellsDefined_ && newPointCount != pointCount_) || (!pointData_.empty() && newPointCount != pointCount_) || (!cellsDefined_ && !cellData_.empty() && newPointCount != cellCount_))
    {
        throw std::logic_error("Cannot change the VTU point count after adding associated data.");
    }

    const std::vector<Real> flattened = flattenVectors(points);
    points_ = makeDataArray(vtkType<Real>(), "", 3, flattened.data(), flattened.size() * sizeof(Real));
    pointCount_ = newPointCount;
    if (!cellsDefined_)
    {
        cellCount_ = pointCount_;
    }
}

void vtuWriter::setCells(const std::vector<int>& connectivity, const std::vector<int>& offsets, const std::vector<unsigned char>& types)
{
    static_assert(sizeof(int) == 4, "VTU Int32 cells require a 32-bit int type.");
    if (offsets.size() != types.size() || (!offsets.empty() && offsets.back() != static_cast<int>(connectivity.size())) || (offsets.empty() && !connectivity.empty()))
    {
        throw std::invalid_argument("Invalid VTU cell storage sizes.");
    }

    int previousOffset = 0;
    for (int offset : offsets)
    {
        if (offset <= previousOffset || offset > static_cast<int>(connectivity.size()))
        {
            throw std::invalid_argument("VTU cell offsets must be strictly increasing.");
        }
        previousOffset = offset;
    }
    for (int pointIndex : connectivity)
    {
        if (pointIndex < 0 || pointIndex >= pointCount_)
        {
            throw std::out_of_range("VTU cell connectivity contains an invalid point index.");
        }
    }
    for (unsigned char type : types)
    {
        if (type == 0)
        {
            throw std::invalid_argument("VTU cell type must not be zero.");
        }
    }

    const int newCellCount = static_cast<int>(types.size());
    if (!cellData_.empty() && newCellCount != cellCount_)
    {
        throw std::logic_error("Cannot change the VTU cell count after adding cell data.");
    }
    connectivity_ = makeDataArray("Int32", "connectivity", 1, connectivity.data(), connectivity.size() * sizeof(int));
    offsets_ = makeDataArray("Int32", "offsets", 1, offsets.data(), offsets.size() * sizeof(int));
    types_ = makeDataArray("UInt8", "types", 1, types.data(), types.size() * sizeof(unsigned char));
    cellCount_ = newCellCount;
    cellsDefined_ = true;
}

void vtuWriter::setVertexCells()
{
    if (!cellData_.empty() && cellCount_ != pointCount_)
    {
        throw std::logic_error("Vertex cells require one cell-data tuple per point.");
    }
    connectivity_ = dataArray{};
    offsets_ = dataArray{};
    types_ = dataArray{};
    cellCount_ = pointCount_;
    cellsDefined_ = false;
}

void vtuWriter::addPointData(const std::string& name, const std::vector<Vec3>& values)
{
    const std::vector<Real> flattened = flattenVectors(values);
    addPointData(name, flattened, 3);
}

void vtuWriter::addCellData(const std::string& name, const std::vector<Vec3>& values)
{
    const std::vector<Real> flattened = flattenVectors(values);
    addCellData(name, flattened, 3);
}

void vtuWriter::write(vtuFormat format) const
{
    if (points_.type_.empty())
    {
        throw std::logic_error("VTU points must be set before writing.");
    }
    createParentDirectory(fileName_);

    dataArray vertexConnectivity;
    dataArray vertexOffsets;
    dataArray vertexTypes;
    const dataArray* outputConnectivity = &connectivity_;
    const dataArray* outputOffsets = &offsets_;
    const dataArray* outputTypes = &types_;
    if (!cellsDefined_)
    {
        std::vector<int> connectivity(pointCount_);
        std::vector<int> offsets(pointCount_);
        std::vector<unsigned char> types(pointCount_, 1);
        for (int pointIndex = 0; pointIndex < pointCount_; ++pointIndex)
        {
            connectivity[pointIndex] = pointIndex;
            offsets[pointIndex] = pointIndex + 1;
        }
        vertexConnectivity = makeDataArray("Int32", "connectivity", 1, connectivity.data(), connectivity.size() * sizeof(int));
        vertexOffsets = makeDataArray("Int32", "offsets", 1, offsets.data(), offsets.size() * sizeof(int));
        vertexTypes = makeDataArray("UInt8", "types", 1, types.data(), types.size() * sizeof(unsigned char));
        outputConnectivity = &vertexConnectivity;
        outputOffsets = &vertexOffsets;
        outputTypes = &vertexTypes;
    }

    std::ofstream output(fileName_, std::ios::binary | std::ios::trunc);
    if (!output)
    {
        throw std::runtime_error("Cannot open VTU file for writing: " + fileName_);
    }

    std::uint64_t appendedOffset = 0;
    std::vector<const dataArray*> appendedArrays;
    const auto writeArray = [&](const dataArray& array)
    {
        writeDataArray(output, array, format, appendedOffset);
        if (format == vtuFormat::binary)
        {
            appendedArrays.push_back(&array);
        }
    };

    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if (format == vtuFormat::binary)
    {
        output << " header_type=\"UInt64\"";
    }
    output << ">\n"
           << "<UnstructuredGrid>\n"
           << "<Piece NumberOfPoints=\"" << pointCount_ << "\" NumberOfCells=\"" << cellCount_ << "\">\n"
           << "<Points>\n";
    writeArray(points_);
    output << "</Points>\n<Cells>\n";
    writeArray(*outputConnectivity);
    writeArray(*outputOffsets);
    writeArray(*outputTypes);
    output << "</Cells>\n";

    if (!pointData_.empty())
    {
        output << "<PointData>\n";
        for (const dataArray& array : pointData_)
        {
            writeArray(array);
        }
        output << "</PointData>\n";
    }
    if (!cellData_.empty())
    {
        output << "<CellData>\n";
        for (const dataArray& array : cellData_)
        {
            writeArray(array);
        }
        output << "</CellData>\n";
    }

    output << "</Piece>\n</UnstructuredGrid>\n";
    if (format == vtuFormat::binary)
    {
        output << "<AppendedData encoding=\"raw\">\n_";
        for (const dataArray* array : appendedArrays)
        {
            writeAppendedData(output, *array);
        }
        output << "\n</AppendedData>\n";
    }
    output << "</VTKFile>\n";
    output.flush();
    if (!output)
    {
        throw std::runtime_error("Failed while writing VTU file: " + fileName_);
    }
}

} // namespace fundem
