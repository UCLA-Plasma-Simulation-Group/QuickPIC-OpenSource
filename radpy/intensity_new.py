import h5py
import time 
import numpy as np
from numpy import cos, linspace, pi, sin
from numpy.linalg import norm
import scipy.constants as sc
import  sys
import os
from mpi4py import MPI 
from scipy import interpolate
#列出计算所需常数 List the constants required for the calculation
c = sc.c#光速 light speed
m = sc.m_e#电子质量 electron mass
e =sc.e #电子电荷 Electron charge
ep = sc.epsilon_0 #epsilon_0
A = e**2/(16*pi**3*c*ep) #积分常数 Integral constant
h = sc.hbar
w0=(4*pi*4e16/(m*1e3))**0.5*e*c*10
def mkdir(path):
    # 引入模块
    # 去除首位空格
    path=path.strip()
    # 去除尾部 \ 符号
    path=path.rstrip("\\")
 
    # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists=os.path.exists(path)
 
    # 判断结果
    if not isExists:
        # 如果不存在则创建目录
        # 创建目录操作函数
        os.makedirs(path) 
 
        print (path+' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print (path+' 目录已存在')
        return False


def pos_x(theta,phi):
    y = np.linspace(theta[0],theta[1],theta[2])
    x = np.linspace(phi[0],phi[1],phi[2])
    posi_x = np.zeros((theta[2],phi[2],3))
    for i in range(phi[2]):
        for j in range(theta[2]):
            posi_x[i,j,:]=[x[i],y[j],(1-x[i]**2-y[j]**2)**0.5]
    return posi_x

def my_interp(z_old,val_old,z_new):
    if val_old.ndim == 1:
        f=interpolate.interp1d(z_old,val_old,kind='linear')
        val_new = f(z_new)
    else:
        val_new = np.zeros((z_new.shape[0],val_old.shape[1]))
        for i in range(val_old.shape[1]):
            f=interpolate.interp1d(z_old,val_old[:,i],kind='linear')
            val_new[:,i] = f(z_new)
    return val_new
def new_orbit(posi,p,t,n,ds):
    vel = p/((p[:,0]**2+p[:,1]**2+p[:,2]**2+1)**0.5).reshape((-1,1))
    p_new= np.array([(p[i,:]+p[i-1,:])/2 for i in range(1,len(t))])
    vel_n = p_new/((p_new[:,0]**2+p_new[:,1]**2+p_new[:,2]**2+1)**0.5).reshape((-1,1))
    delta_t = ds/vel_n[:,2].reshape(-1,1)
    acc = np.array([(vel[i,:]-vel[i-1,:])/delta_t[i-1] for i in range(1,len(t))])
    if n !=0:
        #t_0=t-0.5*delta_t
        posi_0 = posi[1:-1,:]- vel[1:-1,:]*0.5*(delta_t[1:])-0.5*acc[1:,:]*(0.5*(delta_t[1:]))**2
        t_0= t[1:-1]-0.5*ds/vel[:,2][1:-1]
        posi_new=posi_0
        vel_new=vel[1:-1,:]
        acc_new=acc[1:,:]
        t_new=t_0.reshape(-1,1)
        dt=delta_t[1:]/n
        for n in range(1,n):
            vel_in=vel[1:-1,:]+acc[1:,:]*(n*dt)
            posi_in=posi_0+0.5*acc[1:,:]*(n*dt)**2+vel[1:-1,:]*(n*dt)
            t_in=t_0.reshape(-1,1)+n*dt
            acc_in=acc[1:,:]
            posi_new=np.append(posi_new,posi_in,axis=0)
            vel_new=np.append(vel_new,vel_in,axis=0)
            acc_new=np.append(acc_new,acc_in,axis=0)
            t_new=np.append(t_new,t_in,axis=0)
        orbit_new = np.hstack((posi_new,vel_new,acc_new))
    else:
        p_new= np.array([(p[i,:]+p[i-1,:])/2 for i in range(1,len(t))])
        vel_n = p_new/((p_new[:,0]**2+p_new[:,1]**2+p_new[:,2]**2+1)**0.5).reshape((-1,1))
        orbit_new = np.hstack((posi[:-1,:],vel_n,acc))        
        t_new=t[:-1]
    return orbit_new,t_new
def new_orbit2(posi,p,t,n,ds):
    vel = p/((p[:,0]**2+p[:,1]**2+p[:,2]**2+1)**0.5).reshape((-1,1))
    p_new= np.array([(p[i,:]+p[i-1,:])/2 for i in range(1,len(t))])
    vel_n = p_new/((p_new[:,0]**2+p_new[:,1]**2+p_new[:,2]**2+1)**0.5).reshape((-1,1))
    delta_t = ds/vel_n[:,2].reshape(-1,1)
    acc = np.array([(vel[i,:]-vel[i-1,:])/delta_t[i-1] for i in range(1,len(t))])
    k=np.array([(acc[i,:]-acc[i-1,:])/delta_t[i-2] for i in range(1,len(acc))])
    if n !=0:
        #t_0=t-0.5*delta_t
        posi_0 = posi[:-1,:]
        t_0= t[:-1]
        posi_new=posi_0
        vel_new=vel_n
        acc_new=acc
        t_new=t_0.reshape(-1,1)
        dt=delta_t[1:]/n
        for n in range(1,n):
            acc_in=acc[:-1,:]+k*(n*dt)
            vel_in=vel_n[:-1,:]+acc[:-1,:]*(n*dt)+0.5*k*(n*dt)**2
            posi_in=posi_0[:-1,:]+0.5*acc[:-1,:]*(n*dt)**2+vel_n[:-1,:]*(n*dt)+1/6*k*(n*dt)**3
            t_in=t_0[:-1].reshape(-1,1)+n*dt
            posi_new=np.append(posi_new,posi_in,axis=0)
            vel_new=np.append(vel_new,vel_in,axis=0)
            acc_new=np.append(acc_new,acc_in,axis=0)
            t_new=np.append(t_new,t_in,axis=0)
        orbit_new = np.hstack((posi_new,vel_new,acc_new))
    else:
        p_new= np.array([(p[i,:]+p[i-1,:])/2 for i in range(1,len(t))])
        vel_n = p_new/((p_new[:,0]**2+p_new[:,1]**2+p_new[:,2]**2+1)**0.5).reshape((-1,1))
        orbit_new = np.hstack((posi[:-1,:],vel_n,acc))        
        t_new=t[:-1]
    return orbit_new,t_new
def radi_change_omega(posi_particle,posi_x,dt):
    posi = posi_particle[:,0:3]
    vel = posi_particle[:,3:6]
    acc = posi_particle[:,6:9]
    intensity = np.zeros((posi_x.shape[0],posi_x.shape[0]))
    for i in range(posi_x.shape[0]):
        for j in range(posi_x.shape[1]):
            vec_n = [posi_x[i,j,:]]
            temp6=np.sum(vec_n*vel,axis=1)
            vec_n=np.repeat(vec_n,repeats = vel.shape[0],axis=0)
            temp3=np.cross(vec_n - vel,acc,axis=1)
            temp4 = np.cross(vec_n,temp3,axis=1)
            temp5 = ((1- temp6)**5).reshape((1,-1))
            temp1 = np.sum(np.sum(temp4**2,axis=1).reshape((1,-1))/temp5)
            intensity[i,j] += temp1*dt
    return intensity
def main(filename):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank==0:
        mkdir('intensity4')
        mkdir('ELOG4')
        theta=[0.001,-0.001,100]
        phi=[-0.001,0.001,100]
        posi_x = pos_x(theta,phi)
        f = h5py.File(filename,"r")
        keys = list(f.keys())
        keys = np.array_split(keys,size)
        omega = np.linspace(1e5,1.35e8,200)
        start = time.process_time()
    else:
        keys=None
    local_keys = comm.scatter(keys, root=0)
    f = h5py.File(filename,"r")
    posi_x = comm.bcast(posi_x if rank == 0 else None, root = 0)
    omega = comm.bcast(omega if rank == 0 else None, root = 0)
    start = comm.bcast(start if rank == 0 else None, root = 0)
    nn = 0
    inten = np.zeros((posi_x.shape[0],posi_x.shape[1]))
    ff=open("./ELOG4/elog-%08i"%rank,"w")
    ff.write("code begin\n")
    ff.close()
    for part in local_keys:
        p1 = f[part]['p1'][()].astype(float)
        p2 = f[part]['p2'][()].astype(float)
        p3 = f[part]['p3'][()].astype(float)
        p = np.vstack((p1,p2,p3))
        p = p.T
        x1 = f[part]['x1'][()]
        x2 = f[part]['x2'][()]
        x3 = f[part]['x3'][()]
        z = f[part]['t'][()]
        t = (x3 + z)
        posi = np.vstack((x1,x2,z)).T
        orbit_new,t_new=new_orbit2(posi,p,t,5,10)
        dt_new = 2
        inten += radi_change_omega(orbit_new,posi_x,dt_new)
        inten_all = comm.reduce(inten, root=0,op=MPI.SUM)
        nn = nn + 1
        if rank==0:
            np.save("./intensity4/inten%ipart"%(nn*size),inten_all)
        ff=open("./ELOG4/elog-%08i"%rank,"a")
        ff.write("part_num=%i \n"%nn)
        ff.close()

            
    
if __name__ == '__main__':
    main('RAW-BEAM-02_2100-tracks.h5')

