# 获取VASP中能带计算数据的python脚本
### 使用前有一些要求：
- 根据体系调整费米能级，一般要改成scf自洽得到的费米能级，而不是能带计算过程中的费米能级
- 根据计算出的能带是否包含了自旋向上和自旋向下两支选择对应的脚本，**而不是根据是否打开了自旋**，因为有时候即使打开了自旋，计算出的能带还是只有一只
- 在脚本的同级目录下必须含有VASP运行产生的**vasprun.xml**文件和计算能带时用的**KPOINTS**文件
### 具体使用说明参考以下文章：[doc](https://www.jun997.xyz/)
