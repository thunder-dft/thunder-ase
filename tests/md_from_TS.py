import random
from random import randint
import numpy as np
import ase.io
from ase import units
from ase.md.nvtberendsen import NVTBerendsen as NVT
from ase.md.velocitydistribution import ZeroRotation
from ase.calculators.socketio import SocketIOCalculator
from thunder_ase.fireball import Fireball
import os
from multiprocessing import Pool

# 总的运行结构数目
total_sample = 1000
nsample = 4
nthread = 4

temperature = 300  # K
dt = 0.25  # fs
md_step = 10000


def reaction_stat(atoms):
    # 判断反应状态，需要自定义
    status = -1  # 反应未终止
    b1 = atoms.get_distance(3, 9)
    b2 = atoms.get_distance(2, 10)
    if b1 > 3.5 and b2 > 3.5:
        status = 0
    if b1 < 1.6 and b2 < 1.6:
        status = 1
    return status


def run_MD(atoms, name):
    """

    :param atoms:
    :param name: log and traj files name
    :return:
    """
    print("Run the dynamics for {} ...".format(name))
    # FIREBALL settings
    Fdata_path = '/home/ren/data/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'  # set Fdata dir
    # Random Sockets
    unixsocket = 'thunder-ase-{:x}'.format(randint(16 ** 3, 16 ** 4 - 1))
    kwargs = {
        'ipi': 1,
        'nstepf': md_step + 1,  # max step
        'inet': 0,
        'host': unixsocket,
        'scf_tolerance_set': 1E-8,
        'iwriteout_cdcoeffs': 1,
        'iwriteout_charges': 1,
    }
    fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
    fireball.prefix = ''
    status = -1  # before run or failed. 0: Reverse, 1: Forward reaction.
    dyn = NVT(atoms,
              temperature_K=temperature,
              taut=10 * dt * units.fs,
              trajectory=name + '.traj',
              logfile=name + '.log',
              timestep=dt * units.fs,
              )
    os.mkdir(name)
    os.chdir(name)
    with SocketIOCalculator(fireball, log=None, unixsocket=unixsocket) as calc:
        atoms.calc = calc
        for istep, _ in enumerate(dyn.irun(md_step)):
            # 判断是否到达结束标准，是正向还是反向的产物
            status = reaction_stat(atoms)
            if status != -1:
                break
            # Sets the total angular momentum to zero
            ZeroRotation(atoms, preserve_temperature=True)
    os.chdir('..')
    reaction_status = ("Reverse", "Forward", "Failed")[status]
    print("Finish MD for {} after {} steps with reaction status: {}.".format(name, istep, reaction_status))
    return status


if __name__ == '__main__':
    # 键长列表
    rc_list = np.linspace(3.0, 7, 400)  # reaction coordinate list, here is bondlength list
    rc_idx_list = np.arange(400)
    # 初始键长
    ts_bondlength = 4.6
    rc_idx = rc_idx_list[rc_list > ts_bondlength][0]

    record_dict = {}
    while len(record_dict) <= total_sample:
        rc = rc_list[rc_idx]
        traj_name = '{:0>3d}/md-{:.3f}.xyz'.format(rc_idx, rc)
        traj_pool = ase.io.read('../run_16/'+traj_name, index=':')
        print("{} has {} structures.".format(traj_name, len(traj_pool)))
        name_pool = ['bl-{:.3f}-{:0>3d}'.format(rc, i) for i in range(len(traj_pool))]

        traj_pool_idx = range(len(traj_pool))
        # 随机选择 n = 20 个结构
        samples_idx = random.sample(traj_pool_idx, nsample)
        namelist = [name_pool[i] for i in samples_idx]
        # 排除已经计算过的结构
        repeat_idx = []
        repeat_result = []
        for idx, iname in zip(samples_idx, namelist):
            if iname in record_dict:
                repeat_idx.append(idx)
                repeat_result.append(record_dict[iname])
                print("{} has been calculated, with status: {}".format(
                    iname, ("Reverse", "Forward", "Failed")[record_dict[iname]]))

        samples = [traj_pool[i] for i in samples_idx if i not in repeat_idx]
        namelist = [name_pool[i] for i in samples_idx if i not in repeat_idx]
        # 多线程的队列
        with Pool(processes=nthread) as pool:
            joblist = pool.starmap_async(run_MD, zip(samples, namelist))
            results = joblist.get()

        print(results)
        # 去掉未结束的结果
        namelist = [n for n, r in zip(namelist, results) if r != -1]
        results = [r for r in results if r != -1]
        # 记录已经运行的结构，和结果 > dict
        record_dict.update({namelist[idx]: r for idx, r in enumerate(results)})
        # 计算总的正逆反应的比例
        results = results+repeat_result
        assert len(results) > 0
        rf_ratio = sum([i for i in results]) / len(results)  # ratio of reverse and forward reaction
        print("The ratio of reverse and forward reaction is {}.".format(rf_ratio))
        rc_idx += int((rf_ratio - 0.5) * 10)  # update next rc_idx
        print("Next reaction coordinate is {:0>3d}:{:.3f}".format(rc_idx, rc_list[rc_idx]))
