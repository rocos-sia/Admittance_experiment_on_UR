# 使用
## 1. 静止导纳实验

- **仿真测试**
  ```bash
  #打开rviz
  roslaunch ur_examples control.launch  type:=admittance_Control  is_remoted_control:=false
  ```

  ```bash
  #tcp末端施加[10,20,30,10,10,10]N的力
  rostopic pub /force_torque_control std_msgs/Int64MultiArray "layout:
    dim:
    - label: ''
      size: 0
      stride: 0
    data_offset: 0
  data: [10,20,30,10,10,10]"
  ```


- **实物测试**
  ```bash
  roslaunch ur_examples control.launch  type:=admittance_Control  is_remoted_control:=true
  ```
---
## 2. 运动导纳实验
  ### 2.1 直线-力位混合
  - 仿真
    ```bash
    #打开rviz
    roslaunch ur_examples control.launch  type:=moveL_admittance_Control  is_remoted_control:=false
    ``` 
    ```bash
    #TCP末端施加力
    rostopic pub /force_torque_control std_msgs/Int64MultiArray "layout:
      dim:
      - label: ''
        size: 0
        stride: 0
      data_offset: 0
    data: [0,20,30,10,10,10]"
    ```

  - 实物测试
      ```bash
    roslaunch ur_examples control.launch  type:=moveL_admittance_Control  is_remoted_control:=true
    ``` 