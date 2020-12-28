
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
