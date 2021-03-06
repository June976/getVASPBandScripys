# 获取VASP中能带计算数据的python脚本
### 使用前有一些要求：
- 根据体系调整费米能级，一般要改成scf自洽得到的费米能级，而不是能带计算过程中的费米能级
- 根据计算出的能带是否包含了自旋向上和自旋向下两支选择对应的脚本，**而不是根据是否打开了自旋**，因为有时候即使打开了自旋，计算出的能带还是只有一只
- 在脚本的同级目录下必须含有VASP运行产生的**vasprun.xml**文件和计算能带时用的**KPOINTS**文件
### 具体使用说明参考以下文章：[doc](https://www.jun997.xyz/2022/04/18/4ea3d8ef74d8.html)
---
### 更新：添加了两个获取VASP中能带计算数据的bash脚本
### 脚本的正确性验证：[article](https://www.jun997.xyz/2022/04/19/61701185fd17.html)
### 使用脚本注意以下几点：
- bash脚本与python脚本的框架基本是一致的，并且生成的文件两个数据文件格式也是大同小异的
- bash脚本使用的初始输入文件是**OUTCAR**文件和**KPOINTS**文件，这一点与python脚本有区别
