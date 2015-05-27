//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include "device_properties.h"

DeviceProperties::DeviceProperties() {
    m_err = cudaGetDeviceProperties( &m_props, 0 );
}

bool DeviceProperties::good() const {
    return m_err == cudaSuccess;
}

bool DeviceProperties::bad() const {
    return m_err != cudaSuccess;
}

DeviceProperties::~DeviceProperties() {}


std::ostream & operator<<( std::ostream & out, const DeviceProperties & rhs ) {
    if( rhs.bad() ) {
        out << "ERROR\n";
        return out;
    }

    out << "Name: " << rhs.m_props.name << "\n";
    out << "Version: " << rhs.m_props.major << "." << rhs.m_props.minor << "\n";
    out << "computeMode: " << rhs.m_props.computeMode << "\n";
    out << "MaxGridSize: <" << rhs.m_props.maxGridSize[0]
                            << ", " << rhs.m_props.maxGridSize[1]
                            << ", " << rhs.m_props.maxGridSize[2]
                            << ">\n";

    out << "MaxThreadsPerBlock: " << rhs.m_props.maxThreadsPerBlock << "\n";
    out << "sharedMemPerBlock: " << rhs.m_props.sharedMemPerBlock << "\n";

    out << "MaxThreadsDim: <" << rhs.m_props.maxThreadsDim[0]
                                << ", " << rhs.m_props.maxThreadsDim[1] 
                                << ", " << rhs.m_props.maxThreadsDim[2]
                                << ">\n";

    out << "MaxThreadsPerMultiProcessor: " << rhs.m_props.maxThreadsPerMultiProcessor << "\n";
    out << "multiProcessorCount: " << rhs.m_props.multiProcessorCount << "\n";
    out << "MemoryBusWidth: " << rhs.m_props.memoryBusWidth << "\n";
    out << "MemoryClockRate: " << rhs.m_props.memoryClockRate << "\n";

    out << "totalConstMem: " << rhs.m_props.totalConstMem << "\n";
    out << "totalGlobalMem: " << rhs.m_props.totalGlobalMem << "\n";
    return out;
}
