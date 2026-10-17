/**
 * @file HostContainerBindingTest.cpp
 * @brief Checks stable host membership after incremental and bulk vector changes.
 */
#include "data/HostAoSDeviceSoA.h"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{

struct TestObject : fundem::containerObject {
    int value{0};
    struct ValueField;
    using DeviceLayout = fundem::deviceLayout<ValueField>;
};

struct TestObject::ValueField {
    static constexpr auto member = &TestObject::value;
};

using Storage = fundem::hostAoSDeviceSoA<TestObject>;

void require(bool condition, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void checkBindings(const Storage& storage)
{
    const auto& objects = storage.host();
    for (std::size_t index = 0; index < objects.size(); ++index)
    {
        require(objects[index].containerIndex() == static_cast<int>(index), "Host objects must retain their exact container indices.");
    }
}

void testAppendAndRelocation()
{
    Storage growing;
    for (int index = 0; index < 256; ++index)
    {
        growing.host().emplace_back();
        checkBindings(growing);
    }

    Storage storage;
    storage.host().reserve(128);
    for (int index = 0; index < 64; ++index)
    {
        storage.host().emplace_back();
        checkBindings(storage);
    }

    auto& objects = storage.host();
    objects.resize(100);
    checkBindings(storage);
    objects.reserve(1024);
    require(!objects.front().isInContainer(), "Relocated objects must start without copied membership.");
    checkBindings(storage);
    objects.resize(20);
    checkBindings(storage);
    objects.resize(80);
    checkBindings(storage);

    const TestObject copied = storage.host()[10];
    require(!copied.isInContainer(), "An object copy must not inherit container membership.");
    const TestObject moved = std::move(storage.host()[10]);
    require(!moved.isInContainer(), "An object move must not inherit container membership.");
    checkBindings(storage);

    storage.host()[20] = copied;
    require(storage.host()[20].containerIndex() == 20, "Assignment must preserve destination membership.");
    objects.clear();
    objects.assign(40, copied);
    checkBindings(storage);

    Storage assigned;
    assigned.host() = objects;
    checkBindings(assigned);
    assigned.host().assign(12, copied);
    checkBindings(assigned);
    assigned.host().assign(60, copied);
    checkBindings(assigned);

    Storage movedStorage = std::move(storage);
    checkBindings(movedStorage);
    assigned = std::move(movedStorage);
    checkBindings(assigned);
}

void testConcurrentConstAccess()
{
    Storage storage;
    storage.host().resize(4096);
    const Storage& readOnly = storage;
    std::vector<std::thread> readers;
    for (int threadIndex = 0; threadIndex < 4; ++threadIndex)
    {
        readers.emplace_back([&readOnly] {
            for (int iteration = 0; iteration < 100; ++iteration)
            {
                checkBindings(readOnly);
            }
        });
    }
    for (auto& reader : readers)
    {
        reader.join();
    }
}

void benchmarkAppend(int count)
{
    Storage storage;
    storage.host().reserve(static_cast<std::size_t>(count));
    const auto start = std::chrono::steady_clock::now();
    for (int index = 0; index < count; ++index)
    {
        storage.host().emplace_back();
    }
    checkBindings(storage);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << count << " appends: " << seconds << " s\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--benchmark")
        {
            const int count = std::stoi(argv[2]);
            require(count > 0, "The benchmark count must be positive.");
            benchmarkAppend(count);
            return EXIT_SUCCESS;
        }
        testAppendAndRelocation();
        testConcurrentConstAccess();
        std::cout << "Host container membership checks passed.\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
