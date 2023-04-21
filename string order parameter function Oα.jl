using ITensors
using NDTensors
using ITensors.HDF5
using DelimitedFiles
cd("E:/本科学校文件/毕业设计/毕设数据/关联函数/MPS态")

  N = 60             #系统的格点数
  sites = siteinds("S=1",N)     #生成一组ITensor指标对象，数量为N，特性为自旋1
  K = 1                #哈密顿量中的Kitaev项的系数
  J = 0.25             #哈密顿量中的Heisenberg项的系数
  os = OpSum()         #生成一个用于累计哈密顿量的对象

  for j=1:N-1

    os += J,"Sz",j,"Sz",j+1        #通过循环、求和生成哈密顿量中的Heisenberg项
    os += J,"Sy",j,"Sy",j+1
    os += J,"Sx",j,"Sx",j+1

    if rem(j,2) == 0               #通过循环、求和生成哈密顿量中的Kitaev项，使用if语句进行奇偶判断
        os += K,"Sy",j,"Sy",j+1
      else
        os += K,"Sx",j,"Sx",j+1
      end

  end

  os += J,"Sz",1,"Sz",N            #设置PBC边界条件
  os += J,"Sy",1,"Sy",N
  os += J,"Sx",1,"Sx",N
  os += K,"Sy",1,"Sy",N


  H = MPO(os,sites)                #将哈密顿量变换为MPO态

  psi0 = randomMPS(sites,30)       #随机生成一个初始的MPS态，物理指标为sites，所有几何指标的维数为10

  nsweeps = 30                     #nsweep表示DMRG算法的扫描次数
  maxdim = 500                     #定义每次扫描时，几何指标维数的上限
  cutoff = 1E-12                   #定义每次扫描的截断误差目标

  energy, psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)   #dmrg()函数实现了DMRG算法，返回值为基态能量和波函数
  zzcorrx = correlation_matrix(psi,"Sx","Sx")           #计算关联函数

  psix = psi                                           #psix为计算Oα算符期望值过程中的右矢，初始状态时，psix等于系统基态psi
  sx50 = OpSum()                                  
  sx50 += "Sx",50                                      #定义Oα算符中的自旋算符Sx50，由于不同格点的自旋算符对易，自旋算符作用到MPS态上的顺序并不会影响计算结果
  SX50 = MPO(sx50,sites)                               #将自旋算符Sx50变换为MPO态
  psix = apply(SX50,psix)                              #将Sx50自旋算符作用到系统的态函数psix上
  psix = truncate!(psix; cutoff, maxdim)               #对apply函数计算得到的psix态进行剪裁，防止bond dimension过大

  
  i = 49
  while i > 1                                          #使用循环将第2个到第49个节点上的算符exp(im*pi*Sxi)作用到系统的态函数psix上，并使用同样的方法进行剪裁
    sxk = OpSum()
    sxk += "Id",i
    sxk -= 2,"Sx2",i
    SXK = MPO(sxk,sites)
    psix = apply(SXK,psix)
    psix = truncate!(psix; cutoff, maxdim)
    i -= 1
  end


  sx1 = OpSum()                                      #将第一个节点上的自旋算符Sx1作用到系统的态函数psix上，并使用同样的方法进行剪裁
  sx1 += "Sx",1
  SX1 = MPO(sx1,sites)
  psix = apply(SX1,psix)
  psix = truncate!(psix; cutoff, maxdim)

  OalphaX = real(inner(psi,psix))                   #通过计算Oα算符作用后的态函数psix与系统基态psi的内积，得到Oα算符的期望值

