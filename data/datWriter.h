/**
 * @file datWriter.h
 * @brief Declares a column-oriented text writer for simulation histories.
 */
#pragma once

#include "math/Constants.h"

#include <fstream>
#include <initializer_list>
#include <string>
#include <vector>

namespace fundem
{

/** Appends fixed-width numeric rows to a named-column DAT file. */
class datWriter
{
public:
    using Real = math::Real;

    /**
     * Opens a DAT file and writes its column header when required.
     * @param fileName Output path.
     * @param columnNames Ordered names defining the row width.
     * @param appendExisting Preserve and append to an existing file when true.
     * @throws std::runtime_error when the file cannot be opened.
     */
    datWriter(std::string fileName, std::vector<std::string> columnNames, bool appendExisting = true);
    ~datWriter() = default;
    datWriter(const datWriter&) = delete;
    datWriter& operator=(const datWriter&) = delete;
    datWriter(datWriter&&) noexcept = default;
    datWriter& operator=(datWriter&&) noexcept = default;

    const std::string& fileName() const noexcept { return fileName_; }
    int columnCount() const noexcept { return static_cast<int>(columnNames_.size()); }

    /**
     * Appends one numeric row.
     * @param values Values in the same order as the declared columns.
     * @throws std::invalid_argument when the row width differs from `columnCount()`.
     */
    void appendRow(const std::vector<Real>& values);
    /** Convenience overload of `appendRow` for an initializer list. */
    void appendRow(std::initializer_list<Real> values);
    /** Flushes buffered text to the underlying file. */
    void flush();

private:
    /** Emits the fixed column-name header for a newly created DAT file. */
    void writeHeader();

    std::string fileName_;                 ///< Output path.
    std::vector<std::string> columnNames_; ///< Ordered schema written in the header.
    std::ofstream output_;                 ///< Owned output stream.
};

} // namespace fundem
