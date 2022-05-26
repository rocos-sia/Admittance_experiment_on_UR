# 使用
## 1. 导纳实验

### 1.1**仿真测试**
打开rviz
```bash
roslaunch ur_examples control.launch  type:=admittance_Control  real_robot:=false
```
tcp末端施加[10,20,30,10,10,10]N的力
```bash
rostopic pub /force_torque_control std_msgs/Int64MultiArray "layout:
  dim:
  - label: ''
    size: 0
    stride: 0
  data_offset: 0
data: [10,20,30,10,10,10]"
```


### 1.2**实物测试**
```bash
roslaunch ur_examples control.launch  type:=admittance_Control  real_robot:=true
```
