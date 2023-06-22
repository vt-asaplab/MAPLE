# MAPLE: A Metadata-Hiding Policy-Controllable Encrypted Search Platform with Minimal Trust

![x86](https://github.com/vt-asaplab/MAPLE/blob/main/OMAT/EMP_Toolkit/emp-tool/utils/workflows/x86/badge.svg)
![arm](https://github.com/vt-asaplab/MAPLE/blob/main/OMAT/EMP_Toolkit/emp-tool/utils/workflows/arm/badge.svg)

This is our full implementation for our [MAPLE paper].

**WARNING**: This is an academic proof-of-concept prototype and has not received a careful code review. This implementation is NOT ready for production use.

# Required Libraries

1. [ZeroMQ](https://github.com/zeromq/cppzmq/releases/tag/v4.8.1)

2. [EMP-Toolkit](https://github.com/emp-toolkit/emp-agmpc)

You can run the script file **auto_setup.sh** to automatically install the required libraries and build source code. 
```
sudo ./auto_setup.sh
```

## Testing

1. Launch server:
```
$ cd Server
$ ./Server <Server_ID> <Server_Port> 
```

For example, we launch server 1:
```
./Server 1 12345
```
Then launch server 2:
```
./Server 2 12345
```

2. Launch client:
```
$ cd Client
$ ./Client [-b <Bloom_Filter_Size>] [-n <Number_of_Documents>]
```

The default parameters: Bloom_Filter_Size is 16384, Number_of_Documents is 65536. 

For example: 
```
./Client -b 1120 -n 1024
```
