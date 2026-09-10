/**
 * @file vtuWriter.h
 * @brief Declares the low-level ASCII and appended-binary VTU writer.
 */
#pragma once

#include "math/Vector3.h"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <type_traits>
#include <vector>

namespace fundem
{

/** Serialization format for VTU numeric arrays. */
enum class vtuFormat
{
    /** Raw appended binary data. */
    binary,
    ascii
};

/** Builds one VTK XML unstructured-grid file in memory and writes it to disk. */
class vtuWriter
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    /**
     * Creates an empty writer targeting @p fileName.
     * @param fileName Destination VTU path.
     */
    explicit vtuWriter(std::string fileName);

    const std::string& fileName() const noexcept { return fileName_; }
    int pointCount() const noexcept { return pointCount_; }
    int cellCount() const noexcept { return cellCount_; }

    /** Defines point coordinates and establishes the required point-data length. */
    void setPoints(const std::vector<Vec3>& points);
    /** Defines arbitrary VTK cell connectivity, offsets, and cell types. */
    void setCells(const std::vector<int>& connectivity, const std::vector<int>& offsets, const std::vector<unsigned char>& types);
    /** Creates one VTK vertex cell for every point. */
    void setVertexCells();

    /**
     * Adds a scalar or packed-component array associated with points.
     * @tparam Value Arithmetic storage type.
     * @param name VTK array name.
     * @param values Flat values array.
     * @param componentCount Components belonging to each point.
     */
    template <class Value> void addPointData(const std::string& name, const std::vector<Value>& values, int componentCount = 1)
    {
        appendDataArray(pointData_, vtkType<Value>(), name, componentCount, values.data(), values.size() * sizeof(Value), static_cast<int>(values.size()), pointCount_);
    }

    /** Adds a three-component vector array associated with points. */
    void addPointData(const std::string& name, const std::vector<Vec3>& values);

    /**
     * Adds a scalar or packed-component array associated with cells.
     * @tparam Value Arithmetic storage type.
     * @param name VTK array name.
     * @param values Flat values array.
     * @param componentCount Components belonging to each cell.
     */
    template <class Value> void addCellData(const std::string& name, const std::vector<Value>& values, int componentCount = 1)
    {
        appendDataArray(cellData_, vtkType<Value>(), name, componentCount, values.data(), values.size() * sizeof(Value), static_cast<int>(values.size()), cellCount_);
    }

    /** Adds a three-component vector array associated with cells. */
    void addCellData(const std::string& name, const std::vector<Vec3>& values);
    /** Writes raw appended binary by default; pass vtuFormat::ascii for readable XML values. */
    void write(vtuFormat format = vtuFormat::binary) const;

private:
    /** Type-erased numeric field stored until serialization. */
    struct dataArray {
        std::string type_;                ///< VTK scalar type name.
        std::string name_;                ///< User-visible array name.
        int componentCount_{1};           ///< Values stored per entity.
        std::vector<std::uint8_t> bytes_; ///< Native binary representation of all values.
    };

    /** Maps one supported arithmetic C++ type to its VTK XML scalar name. */
    template <class Value> static const char* vtkType()
    {
        using Type = std::remove_cv_t<Value>;
        static_assert(std::is_arithmetic_v<Type> && !std::is_same_v<Type, bool>, "VTU arrays require a non-boolean arithmetic value type.");

        if constexpr (std::is_same_v<Type, float>)
        {
            return "Float32";
        }
        else if constexpr (std::is_same_v<Type, double>)
        {
            return "Float64";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_signed_v<Type> && sizeof(Type) == 1)
        {
            return "Int8";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_unsigned_v<Type> && sizeof(Type) == 1)
        {
            return "UInt8";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_signed_v<Type> && sizeof(Type) == 2)
        {
            return "Int16";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_unsigned_v<Type> && sizeof(Type) == 2)
        {
            return "UInt16";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_signed_v<Type> && sizeof(Type) == 4)
        {
            return "Int32";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_unsigned_v<Type> && sizeof(Type) == 4)
        {
            return "UInt32";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_signed_v<Type> && sizeof(Type) == 8)
        {
            return "Int64";
        }
        else if constexpr (std::is_integral_v<Type> && std::is_unsigned_v<Type> && sizeof(Type) == 8)
        {
            return "UInt64";
        }
        else
        {
            static_assert(sizeof(Type) == 0, "Unsupported VTU value type.");
        }
    }

    /** Copies one typed byte range into type-erased deferred storage. */
    static dataArray makeDataArray(const char* type, const std::string& name, int componentCount, const void* data, std::size_t byteCount);
    /** Writes an array declaration or inline ASCII values and advances binary offset. */
    static void writeDataArray(std::ostream& output, const dataArray& array, vtuFormat format, std::uint64_t& appendedOffset);
    /** Decodes native numeric bytes and emits XML text values. */
    static void writeASCIIData(std::ostream& output, const dataArray& array);
    /** Writes one length-prefixed raw array to the appended-data section. */
    static void writeAppendedData(std::ostream& output, const dataArray& array);
    /** Validates tuple count and appends one deferred point/cell field. */
    static void appendDataArray(std::vector<dataArray>& arrays,
                                const char* type,
                                const std::string& name,
                                int componentCount,
                                const void* data,
                                std::size_t byteCount,
                                int valueCount,
                                int entityCount);

    std::string fileName_;             ///< Destination path.
    int pointCount_{0};                ///< Number of mesh points.
    int cellCount_{0};                 ///< Number of mesh cells.
    bool cellsDefined_{false};         ///< Whether connectivity has been supplied.
    dataArray points_;                 ///< Point-coordinate array.
    dataArray connectivity_;           ///< Flattened cell connectivity.
    dataArray offsets_;                ///< Cumulative connectivity offsets.
    dataArray types_;                  ///< VTK cell-type identifiers.
    std::vector<dataArray> pointData_; ///< User-selected point fields.
    std::vector<dataArray> cellData_;  ///< User-selected cell fields.
};

} // namespace fundem
