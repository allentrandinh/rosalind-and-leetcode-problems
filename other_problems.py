
#Fibonacci sequence
#start out with 1 pair of rabbits, each pair takes a month to start having offspring
#input (n,k) with n: number of generations, k: number of pairs of offspring produced each month
def rabbit_count(num_gen,num_off):
    if num_gen == 0:
        return 0
    if num_gen == 1:
        return 1
    count = rabbit_count(num_gen-1,num_off) + num_off*rabbit_count(num_gen-2,num_off)
    return count



