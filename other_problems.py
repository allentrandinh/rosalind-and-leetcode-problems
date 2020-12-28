
#start out with 1 pair of rabbits, each pair takes a month to start having offspring
def fibonacci_rabbits(num_gen,num_off):
    '''
    immortal rabbits
    :param num_gen: number of generation
    :param num_off: number of offspring pairs produced each month by 1 adult rabbit pair.
    :return: number of rabbits
    '''
    if num_gen == 0:
        return 0
    if num_gen == 1:
        return 1
    count = rabbit_count(num_gen-1,num_off) + num_off*rabbit_count(num_gen-2,num_off)
    return count

class MortalFibonacciRabbit:
    '''
    start out with 1 pair of rabbit. Each pair takes a month before they start to produce offsprings.
    rabbits die after certain amount of months
    create an instances by specifying number of generations and survival months
    instance_name.mortal_rabbit() returns number of rabbits at the specified generation
    '''
    def __init__(self,num_generation,survival_month):
        self.num_generation = num_generation
        self.survival_month = survival_month

    @staticmethod
    def before_first_rabbit_die(num_gen, survival_month):
        new_rabbit, old_rabbit = dict[num_gen - 1][1], dict[num_gen - 1][2]
        total_rabbit = new_rabbit + old_rabbit
        return (new_rabbit, old_rabbit, total_rabbit)

    @staticmethod
    def after_first_rabbit_die(num_gen, survival_month):
        new_rabbit, old_rabbit = dict[num_gen - 1][1], dict[num_gen - 1][2] - dict[num_gen - survival_month][0]
        total_rabbit = new_rabbit + old_rabbit
        return (new_rabbit, old_rabbit, total_rabbit)

    @staticmethod
    def rabbit_dictionary_generation(generation, survival_month):
        i = 2
        while i <= survival_month:
            dict.update({i: before_first_rabbit_die(i, survival_month)})
            i += 1
        while i <= generation:
            dict.update({i: after_first_rabbit_die(i, survival_month)})
            i += 1
        return (dict)

    def mortal_rabbit(self):
        dict = {1: (1, 0, 1)}
        result = rabbit_dictionary_generation(generation=self.num_generation, survival_month=self.survival_month)
        return (result[self.num_generation][2])

def dominant_prob(k,m,n):
    '''
    :param k: number homozygous dominant invidividuals in pplt
    :param m: number of heterozygous
    :param n: number of homozygous recessive
    :return: probability that two randomly selected mating organisms will produce an individual possessing a dominant allele
    '''
    pplt_size = k + m + n
    #recessive*recessive=no dominant, recessive*hetero = 1/2, hetero*hetero=1/4
    not_dominant = (n/pplt_size)*((n-1)/(pplt_size-1)) \
                   + (n/pplt_size)*(m/(pplt_size-1)) \
                   + (m/pplt_size)*((m-1)/(pplt_size-1))*1/4
    return(1-not_dominant)


