using ITensors





let
  N = 60       #系统的格点数
  sites = siteinds("S=1",N)    #生成一组ITensor指标对象，数量为N，特性为自旋1
  K = 1               #哈密顿量中的Kitaev项的系数
  J = 0.5             #哈密顿量中的Heisenberg项的系数
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


  os += J,"Sz",1,"Sz",N
  os += J,"Sy",1,"Sy",N
  os += J,"Sx",1,"Sx",N
  os += K,"Sy",1,"Sy",N


  H = MPO(os,sites)              #将哈密顿量变换为MPO态

  psi0 = randomMPS(sites,30)     #随机生成一个初始的MPS态，物理指标为sites，所有几何指标的维数为30

  nsweeps = 30                     #nsweep表示DMRG算法的扫描次数
  maxdim = 500                   #定义每次扫描时，几何指标维数的上限
  cutoff = 1E-16                 #定义每次扫描的截断误差目标

  energy, psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)  #dmrg()函数实现了DMRG算法，返回值为基态能量和波函数

  b = 30
  orthogonalize!(psi, b)      #将psi的正交中心移动到b
  U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))     #施密特分解

  

  print("\n")
  print("\n")
  print("\n")
  print("\n")
  print("\n")
  print("\n")


  for i=1:500        #输出纠缠谱

     print("alpha = ")
     print(i)
     print("      Saa = ")
     print(S[i,i]^2)
     print("\n")

  end

end








  
  