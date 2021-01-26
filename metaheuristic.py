import os
import pickle
import pandas as pd
import numpy as np
import seaborn as sns
import xlwings as xw
import matplotlib.pyplot as plt
import time
import random
sns.set()
sns.set_style('darkgrid')
random.seed(1)

def start_mh(ct_,p_, s_, depots_, back_, pt_, dst_mtx_, indxs_, prb_thrs_ = 0.8):
    pt_ = pt_.sort_values(axis = 0)
    k = 1
    while k>0:
        sln_vctr_ = create_init_sol(depots_, back_, pt_, dst_mtx_, indxs_, prb_thrs_)
        of_, rt_, t_ = compute_of_2(sln_vctr_, dst_mtx_, ct_,p_, s_)
        
        #print(sln_vctr_)
        
        if sum(t_.sort_index(axis = 0) > pt_.sort_index(axis = 0)) == 0 and sum(t_.sort_index(axis = 0) < ct_.sort_index(axis = 0)) == 0:
            k = 0
    
    return sln_vctr_, of_, rt_, t_

def best_insertion(dst_mtx_,depots_, back_, ct_, p_, s_, pt_):
    sln_vctr_ = []
    for i in range(len(depots_)):
        sln_vctr_.append([depots_[i]])
    for idx, i_ in enumerate(sln_vctr_):
        i_.append(back_[idx])
    
    nodes_list_ = list(pt_.index)
    random.shuffle(nodes_list_)
    n_sol_ = sln_vctr_.copy()
    for i_ in nodes_list_:
        n_sol_, n_of_, n_rt_, n_t_ = insert_node(n_sol_, i_, dst_mtx_, ct_, p_, s_, pt_)
        
        
    return n_sol_, n_of_, n_rt_, n_t_

def insert_node(sln_vctr_,n_ins_, dst_mtx_, ct_, p_, s_, pt_):
    
    b_sol_ = []
    b_of_ = np.inf
    b_rt_ = None
    b_t_ = None

    for rdx, route_ in enumerate(sln_vctr_):
        for idx, i_ in enumerate(route_):
            if idx != len(route_)-1:
                n_route_ = route_[:idx+1] + [n_ins_] + route_[idx+1:]

                n_sln_vctr_ = sln_vctr_[:rdx] + [n_route_] + sln_vctr_[rdx+1:]
                #print(n_sln_vctr)
                n_of_, n_rt_, n_t_ = compute_of_2(n_sln_vctr_, dst_mtx_, ct_,p_, s_)

                if sum(n_t_.sort_index(axis = 0) > pt_[n_t_.index].sort_index(axis = 0)) == 0 and sum(n_t_.sort_index(axis = 0) < ct_[n_t_.index].sort_index(axis = 0)) == 0:

                    if n_of_ < b_of_:
                        b_sol_ = n_sln_vctr_.copy()
                        b_of_ = n_of_
                        b_rt_ = n_rt_.copy()
                        b_t_ = n_t_.copy()
    
    return b_sol_, b_of_, b_rt_, b_t_
    

def create_init_sol(depots_, back_, pt_, dst_mtx_, indxs_, prb_thrs_ = 0.8):
    
    if len(depots_) == 1:
        sln_vctr_ = [depots_ + list(pt_.index) + back_]
    
    else:
        
        sln_vctr_ = []
        for i in range(len(depots_)):
            sln_vctr_.append([depots_[i]])
        nodes_from_ = [i[-1] for i in sln_vctr_]
        
        for idx, item in pt_.iteritems():
            rnd = random.random()
            nodo_ = idx
            if rnd <= prb_thrs_:
                
                t_n_ = dst_mtx_.loc[nodes_from_,nodo_].idxmin()
                route_number_ = nodes_from_.index(t_n_)
                sln_vctr_[route_number_].append(nodo_)
                
            else:
                num_nodes_ = [len(i) for i in sln_vctr_]
                min_rt_ = min(num_nodes_)
                min_idx_ = num_nodes_.index(min_rt_)
                sln_vctr_[min_idx_].append(nodo_)
            nodes_from_ = [i[-1] for i in sln_vctr_]
            
        for idx, item in enumerate(sln_vctr_):
            item.append(back_[idx])
        
        
        
    return sln_vctr_



def compute_of_2(sln_vctr_, dst_mtx_, ct_,p_, s_):
    rt_ = pd.Series(dtype='float64')
    t_ = pd.Series(dtype='float64')
    #wt_t_ = pd.Series(dtype = 'float64')
    for route_ in sln_vctr_:
        
        l_ = len(route_)-1
        if l_ > 1:
            for idx_, i_ in enumerate(route_):
                if idx_ == 1:
                    tr_temp_ = dst_mtx_.loc[route_[idx_-1],i_]
                    t_temp_ = ct_.loc[i_] + dst_mtx_.loc[route_[idx_-1],i_]
                    t_.loc[i_] = t_temp_
                    rt_.loc[i_] = tr_temp_
                    #wt_t_.loc[i_] = 0
                elif idx_!= 0 and idx_ != l_:
                    t_temp_ = dst_mtx_.loc[route_[idx_-1],i_] + max(s_.loc[route_[idx_-1]] + t_.loc[route_[idx_-1]], ct_.loc[i_])
                    tr_temp_ = t_temp_ - ct_.loc[i_]
                    rt_.loc[i_] = tr_temp_
                    t_.loc[i_] = t_temp_
    
    of_ = np.sum(rt_*p_)
    
    return of_, rt_, t_



def swap_ns(i_,j_, sln_vctr_):
    
    idx_ = sln_vctr_.index(i_)
    h_ = sln_vctr_[idx_-1]
    k_ = sln_vctr_[idx_+1]
    jdx_ = sln_vctr_.index(j_)
    l_ = sln_vctr_[jdx_-1]
    m_ = sln_vctr_[jdx_+1]    
    
    idx_ = sln_vctr_.index(i_)
    jdx_ = sln_vctr_.index(j_)
    n_solvctr_ = list(sln_vctr_)
    n_solvctr_[idx_] = j_
    n_solvctr_[jdx_] = i_
    
    
    return n_solvctr_


def run_swap_ns(sln_vctr_, obj_fnc_, dst_mtx_, ct_vctr_, p_, pt_, s_):
    #print(time.ctime())
    best_sln_vct_ = sln_vctr_.copy()
    best_obj_fnc_ = obj_fnc_
    _,best_rt_, best_t_ = compute_of_2(sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
    
    for rdx, route in enumerate(sln_vctr_):
        for idx, item in enumerate(route):
            for jdx, jtem in enumerate(route[idx+1:-1]):
                
                if idx != 0:
                    #print(item, route)
                    n_route_ = swap_ns(item,jtem, route)
                    n_sln_vctr_ = sln_vctr_[:rdx] + [n_route_] + sln_vctr_[rdx+1:]
                    
                    n_obj_fnc_, n_rt_, n_t_ = compute_of_2(n_sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
                    
                    
                    #n_obj_fnc_, n_at_,n_t_ = compute_cumulative_of(dst_mtx_, n_sln_mtx_, wt_mtx_, ct_vctr_, n_ord_vct_, p_)
                    #print(n_obj_fnc_)
                    if sum(n_t_.loc[pt_.index] > pt_) == 0 and sum(n_t_.loc[ct_vctr_.index] < ct_vctr_) == 0:
                        if n_obj_fnc_ < best_obj_fnc_:
                            best_sln_vct_ = n_sln_vctr_.copy()
                            best_obj_fnc_ = n_obj_fnc_
                            best_rt_ = n_rt_
                            best_t_ = n_t_
    #print(time.ctime())
    return best_sln_vct_, best_obj_fnc_, best_rt_, best_t_

def relocate_ns(i_, j_, sln_vctr_):
    
    idx_ = sln_vctr_.index(i_)
    jdx_ = sln_vctr_.index(j_)
    hdx_ = idx_-1
    h_ = sln_vctr_[hdx_]
    
    
    head_ = list(sln_vctr_[:idx_])
    mid_ = list(sln_vctr_[idx_+1:jdx_+1])
    tail_ = list(sln_vctr_[jdx_+1:])
    n_solvctr_ = head_+mid_+[i_]+tail_
    return n_solvctr_



def run_relocate_ns(sln_vctr_, obj_fnc_, dst_mtx_, ct_vctr_, p_, pt_, s_):
    #print(time.ctime())
    best_sln_vct_ = sln_vctr_.copy()
    best_obj_fnc_ = obj_fnc_
    _,best_rt_, best_t_ = compute_of_2(sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
    
    for rdx, route in enumerate(sln_vctr_):
        for idx, item in enumerate(route):
            for jdx, jtem in enumerate(route[idx+1:-1]):
                
                if idx != 0:
                    #print(item, route)
                    n_route_ = relocate_ns(item,jtem, route)
                    n_sln_vctr_ = sln_vctr_[:rdx] + [n_route_] + sln_vctr_[rdx+1:]
                    
                    n_obj_fnc_, n_rt_, n_t_ = compute_of_2(n_sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
                    
                    
                    #n_obj_fnc_, n_at_,n_t_ = compute_cumulative_of(dst_mtx_, n_sln_mtx_, wt_mtx_, ct_vctr_, n_ord_vct_, p_)
                    #print(n_obj_fnc_)
                    if sum(n_t_.loc[pt_.index] > pt_) == 0 and sum(n_t_.loc[ct_vctr_.index] < ct_vctr_) == 0:
                        if n_obj_fnc_ < best_obj_fnc_:
                            best_sln_vct_ = n_sln_vctr_.copy()
                            best_obj_fnc_ = n_obj_fnc_
                            best_rt_ = n_rt_
                            best_t_ = n_t_
    #print(time.ctime())
    return best_sln_vct_, best_obj_fnc_, best_rt_, best_t_

def ir_swap_ns(ri_, rj_, i_,j_, sln_vctr_):
    
    idx_ = sln_vctr_[ri_].index(i_)
    h_ = sln_vctr_[ri_][idx_-1]
    k_ = sln_vctr_[ri_][idx_+1]
    
    jdx_ = sln_vctr_[rj_].index(j_)
    l_ = sln_vctr_[rj_][jdx_-1]
    m_ = sln_vctr_[rj_][jdx_+1]
    
    
    idx_ = sln_vctr_[ri_].index(i_)
    jdx_ = sln_vctr_[rj_].index(j_)
    n_solvctr_ = [[i for i in ii ] for ii in sln_vctr_]
    n_solvctr_[ri_][idx_] = j_
    n_solvctr_[rj_][jdx_] = i_
    
    return n_solvctr_




def run_ir_swap_ns(sln_vctr_, obj_fnc_, dst_mtx_, ct_vctr_, p_, pt_, s_):
    #print(time.ctime())
    best_sln_vct_ = sln_vctr_.copy()
    best_obj_fnc_ = obj_fnc_
    _,best_rt_, best_t_ = compute_of_2(sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
    
    
    for ridx, routei in enumerate(sln_vctr_):
        for rjdx, routej in enumerate(sln_vctr_[ridx+1:]):
    
            for idx, item in enumerate(routei[1:-1]):
                for jdx, jtem in enumerate(routej[1:-1]):

                    if idx != 0:
                        #print(item, route)
                        n_sln_vctr_ = ir_swap_ns(ridx,rjdx+ridx+1,item,jtem, sln_vctr_)

                        n_obj_fnc_, n_rt_, n_t_ = compute_of_2(n_sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)


                        #n_obj_fnc_, n_at_,n_t_ = compute_cumulative_of(dst_mtx_, n_sln_mtx_, wt_mtx_, ct_vctr_, n_ord_vct_, p_)
                        #print(n_obj_fnc_)
                        if sum(n_t_.loc[pt_.index] > pt_) == 0 and sum(n_t_.loc[ct_vctr_.index] < ct_vctr_) == 0:
                            if n_obj_fnc_ < best_obj_fnc_:
                                best_sln_vct_ = n_sln_vctr_.copy()
                                best_obj_fnc_ = n_obj_fnc_
                                best_rt_ = n_rt_
                                best_t_ = n_t_
    #print(time.ctime())
    return best_sln_vct_, best_obj_fnc_, best_rt_, best_t_

def ir_relocate_ns(ri_, rj_, i_, j_, sln_vctr_):
    
    idx_ = sln_vctr_[ri_].index(i_)
    jdx_ = sln_vctr_[rj_].index(j_)
    hdx_ = idx_-1
    h_ = sln_vctr_[ri_][hdx_]
    
    n_solvctr_ = [[i for i in ii ] for ii in sln_vctr_]
    
    n_solvctr_[ri_] = n_solvctr_[ri_][:idx_] + sln_vctr_[ri_][idx_+1:]
    n_solvctr_[rj_] = n_solvctr_[rj_][:jdx_+1] + [i_] + sln_vctr_[rj_][jdx_+1:]
    
    
    return n_solvctr_

def run_ir_relocate_ns(sln_vctr_, obj_fnc_, dst_mtx_, ct_vctr_, p_, pt_, s_):
    #print(time.ctime())
    best_sln_vct_ = sln_vctr_.copy()
    best_obj_fnc_ = obj_fnc_
    _,best_rt_, best_t_ = compute_of_2(sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
    
    
    for ridx, routei in enumerate(sln_vctr_):
        for rjdx, routej in enumerate(sln_vctr_):
            
            if ridx != rjdx:
                
                for idx, item in enumerate(routei[1:-1]):
                    for jdx, jtem in enumerate(routej[:-1]):

                        if idx != 0:
                            #print(item, route)
                            n_sln_vctr_ = ir_relocate_ns(ridx,rjdx,item,jtem, sln_vctr_)

                            n_obj_fnc_, n_rt_, n_t_ = compute_of_2(n_sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)


                            #n_obj_fnc_, n_at_,n_t_ = compute_cumulative_of(dst_mtx_, n_sln_mtx_, wt_mtx_, ct_vctr_, n_ord_vct_, p_)
                            #print(n_obj_fnc_)
                            if sum(n_t_.loc[pt_.index] > pt_) == 0 and sum(n_t_.loc[ct_vctr_.index] < ct_vctr_) == 0:
                                if n_obj_fnc_ < best_obj_fnc_:
                                    best_sln_vct_ = n_sln_vctr_.copy()
                                    best_obj_fnc_ = n_obj_fnc_
                                    best_rt_ = n_rt_
                                    best_t_ = n_t_
    #print(time.ctime())
    return best_sln_vct_, best_obj_fnc_, best_rt_, best_t_


def VND_heur(sln_vctr_, obj_fnc_, dst_mtx_, ct_vctr_, p_, pt_, s_, heuristics_list = ['ir-relocate','ir-swap','relocate', 'swap']):
    #print(time.ctime())
    k_ = len(heuristics_list)
    i_ = 0
    
    best_sln_vct_ = sln_vctr_
    best_of_ = obj_fnc_
    _,best_rt_, best_t_ = compute_of_2(sln_vctr_, dst_mtx_, ct_vctr_ ,p_, s_)
    
    ii_ = 0
    while i_ <= k_ - 1:
        ii_ = ii_+1
        
        if heuristics_list[i_] == 'ir-relocate' and len(sln_vctr_) >1:
            n_sln_vct_, n_of_, n_rt_, n_t_ = run_ir_relocate_ns(best_sln_vct_, best_of_, dst_mtx_, ct_vctr_, p_, pt_, s_)
            
        elif heuristics_list[i_] == 'ir-swap' and len(sln_vctr_) > 1:
            n_sln_vct_, n_of_, n_rt_, n_t_ = run_ir_swap_ns(best_sln_vct_, best_of_, dst_mtx_, ct_vctr_, p_, pt_, s_)
            
        elif heuristics_list[i_] == 'relocate':
            n_sln_vct_, n_of_, n_rt_, n_t_ = run_relocate_ns(best_sln_vct_, best_of_, dst_mtx_, ct_vctr_, p_, pt_, s_)
        
        elif heuristics_list[i_] == 'swap':
            n_sln_vct_, n_of_, n_rt_, n_t_ = run_swap_ns(best_sln_vct_, best_of_, dst_mtx_, ct_vctr_, p_, pt_, s_)
              
        if n_of_ < best_of_:
            best_sln_vct_ = n_sln_vct_.copy()
            best_of_ = n_of_
            best_rt_ = n_rt_
            best_t_ = n_t_
            #print(heuristics_list[i_], n_of_)
    
            i_ = 0
        else:
            i_ = i_ + 1
        #print(ii_, best_of_, heuristics_list[i_])
    #print(time.ctime())        
    return best_sln_vct_,  best_of_, best_rt_, best_t_


