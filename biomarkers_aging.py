import numpy as np
import random
import math
import sys
import scipy
from scipy.optimize import minimize
import sys
#replicate = sys.argv[1]

num_loci = 20000#the total number of loci in the sumulation
num_samples_per_catagory = 1 #the number of individuals per target cohort group
#p_locus_age = np.append(np.ones((num_loci/2)) * 0.001, np.ones((num_loci/2)) * 0.0005)

def sim_ind_binom(p_locus_age, locus_effect_size, ages_of_interest, num_loci, frac__mort_loci_dependent):
    initial_mortality = math.exp(-9.575)
    dead = 0
    time = 1
    locus_status = np.zeros((num_loci))
    locus_aged_date = np.zeros((num_loci))
    sum_loci_aged = 0
    sum_age_effect = 0
    mortality = initial_mortality
    mort_at_ages_of_interest = np.zeros((len(ages_of_interest)+1))
    record_num = 0
    total_effect_size = sum(locus_effect_size)
    mean_effect_size = total_effect_size / float(len(locus_effect_size))
    mean_p_deg = sum(p_locus_age) / float(len(p_locus_age))
    exponent = (0.0807 / math.log(1-1/76))
    sum_loci_eff = np.zeros((120))
    mortality = np.zeros((len(sum_loci_eff)))
    locus = 0
    locus_aged_date = np.random.negative_binomial(1,p_locus_age, len(p_locus_age))
    if frac__mort_loci_dependent >0:
        for deg_date in locus_aged_date:
            if deg_date < len(sum_loci_eff) and locus_effect_size[locus] > 0:
                affected_loci = sum_loci_eff[deg_date:len(sum_loci_eff)] + locus_effect_size[locus] * np.ones((len(sum_loci_eff)-deg_date))
                sum_loci_eff[deg_date:len(sum_loci_eff)] = affected_loci
            locus+=1
    for age in range(1,len(mortality)):
        gompertz_mort = math.exp(0.0807*age-9.575)
        #print(sum_loci_eff[age])
        loci_induced_mortality = initial_mortality *(1 - sum_loci_eff[age]/ total_effect_size)**exponent
        mortality[age] = (1-frac__mort_loci_dependent) * gompertz_mort + frac__mort_loci_dependent * loci_induced_mortality
    age = 0
    while np.random.random() > mortality[age] and age < len(sum_loci_eff)-1:
        if record_num < len(ages_of_interest) and age == ages_of_interest[record_num]:
            mort_at_ages_of_interest[record_num] = mortality[age]
            record_num += 1
        age+=1
        #print (time, sum_age_effect, mortality)
    return(age, locus_aged_date, mort_at_ages_of_interest)


def simulate_cohorts(p_locus_age, locus_effect_size, age_categories, frac__mort_loci_dependent, replicate):
    mean_age = 0
    rep = str(replicate)
    mean_lifepsan_remaining = 0
    mean_loci_status = np.zeros((num_loci))
    loci_statusXage = np.zeros((num_loci))
    loci_statusXlifepsan_remaining = np.zeros((num_loci))
    ageXage = 0
    lifepsan_remainingXlifepsan_remaining = 0
    loci_statusXloci_status = np.zeros((num_loci))
    loci_statusXmort_rate = np.zeros((num_loci))
    mean_mort_rate = 0
    print_all_ind = "yes"
    if print_all_ind == "yes":
        filename = "loci_status" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age))) + "_" + rep + ".txt"
        loci_status = open(filename , "w")
        filename = "ages" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age)))+  "_" + rep + ".txt"
        ages = open(filename , "w")
        filename = "lifespan_remaining" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age)))+  "_" + rep +".txt"
        lifespan_remaining_file = open(filename , "w")
        filename = "aging_rate" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age)))+  "_" + rep +".txt"
        aging_rate_file = open(filename , "w")

        loci_info = open("loci_info" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age))) +  "_" + rep +".txt" , "w")
        locus = 0
        for p_deg in p_locus_age:
            print((1-p_deg) / p_deg, locus_effect_size[locus], sep = "\t", file = loci_info)
            locus +=1
        loci_info.close()

    #    values = np.zeros((num_samples_per_catagory*len(age_categories), num_loci))
    total_samples = 0
    print_loci_stats = 0
    for age_sampled in age_categories:
        print (age_sampled)
        successful_samples = 0
        while successful_samples < num_samples_per_catagory:
            age_rate_type = 0.0807#np.random.normal(0.0807, 0.0807*0.1)
            for locus in range(0, num_loci):
                p_locus_age[locus] = 1 - math.exp(age_rate_type * math.log(1 - 1/76) / 0.0807)
            simulation = sim_ind_binom(p_locus_age, locus_effect_size, [age_sampled], num_loci, frac__mort_loci_dependent)
            age_at_death = simulation[0]
            #print(age_at_death)
            ages_of_degredation = simulation[1]
            #print(successful_samples, age_at_death)
            mort_rate = simulation[2]
            if age_at_death > age_sampled and  age_at_death - age_sampled < 26:
                mean_age += age_sampled
                ageXage += age_sampled**2
                mean_mort_rate += mort_rate[0]
                mean_lifepsan_remaining += (age_at_death - age_sampled)
                lifepsan_remainingXlifepsan_remaining += (age_at_death - age_sampled)**2
                locus = 0
                num_degraded_loci = 0
                degraded_loci_eff = 0
                for date in ages_of_degredation:
                    if date < age_sampled and date > 0:
                        degraded_loci_eff += locus_effect_size[locus]
                        if locus_effect_size[locus] > 0:
                            num_degraded_loci += 1
                        if print_all_ind == "yes":
                            print(1, end = "\t", file = loci_status)
                    elif print_all_ind == "yes":
                        print(0, end = "\t", file = loci_status)

                        if print_loci_stats == 1:
                            mean_loci_status[locus]+= 1
                            loci_statusXloci_status[locus] += 1
                            loci_statusXage[locus] += age_sampled
                            loci_statusXlifepsan_remaining[locus] = (age_at_death - age_sampled)
                            loci_statusXmort_rate[locus] += mort_rate[0]
                    locus += 1

                print("0", end = "\n", file = loci_status)

                if print_all_ind == "yes":
                    print(age_sampled, file = ages)
                    print(age_rate_type, file = aging_rate_file)
                    lifespan_remaining = age_at_death - age_sampled
                    print(lifespan_remaining, file = lifespan_remaining_file)
                    #print(lifespan_remaining)

                successful_samples +=1
                total_samples += 1
    if(print_loci_stats == 1):
        cov_loci_age = 0
        cov_loci_lifespan_remaining = 0
        cov_loci_mort = 0
        var_loci = 0
        mean_age = mean_age / (num_samples_per_catagory*len(age_categories))
        mean_lifepsan_remaining = mean_lifepsan_remaining / (num_samples_per_catagory*len(age_categories))
        mean_mort_rate = mean_mort_rate/ (num_samples_per_catagory*len(age_categories))
        var_age = ageXage / (num_samples_per_catagory*len(age_categories)) - mean_age**2
        var_lifespan_remaining = lifepsan_remainingXlifepsan_remaining / (num_samples_per_catagory*len(age_categories)) - mean_lifepsan_remaining**2

        filename = "loci_aging" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age)))+ ".txt"
        f = open(filename , "w")
        print("deg_freq", "locus_effect_size","mean_loci_status","var_loci", "cov_loci_mort","cov_loci_age","corr_loci_age","cov_loci_lifespan_remaining","corr_loci_lifespan_remaining", sep = "\t", file = f)
        for locus in range(0,num_loci):
            mean_loci_status[locus] = mean_loci_status[locus] / (num_samples_per_catagory*len(age_categories))
            var_loci = loci_statusXloci_status[locus] / (num_samples_per_catagory*len(age_categories)) - mean_loci_status[locus]**2
            cov_loci_age = loci_statusXage[locus] / (num_samples_per_catagory*len(age_categories)) - mean_loci_status[locus] * mean_age
            cov_loci_lifespan_remaining = loci_statusXlifepsan_remaining[locus] / (num_samples_per_catagory*len(age_categories)) - mean_loci_status[locus] * mean_lifepsan_remaining
            cov_loci_mort = loci_statusXmort_rate[locus] / (num_samples_per_catagory*len(age_categories)) - mean_mort_rate * mean_loci_status[locus]
            corr_loci_age = 0
            corr_loci_lifespan_remaining = 0
            if var_loci > 0:
                corr_loci_age = cov_loci_age / (var_age**0.5 * var_loci**0.5)
                corr_loci_lifespan_remaining = cov_loci_lifespan_remaining / (var_lifespan_remaining**0.5 * var_loci**0.5)

            print(1 / p_locus_age[locus], locus_effect_size[locus], mean_loci_status[locus], var_loci, cov_loci_mort, cov_loci_age, corr_loci_age, cov_loci_lifespan_remaining, corr_loci_lifespan_remaining, sep = "\t", file = f)
        f.close()
    ages.close()
    loci_status.close()
        #print(mean_lifespan, locus_effect_size[locus], num_degraded[locus], num_not_degraded[locus], mean_date_aged[locus]/num_degraded[locus], sep = "\t", file = f)

def simulate_pop(p_locus_age, locus_effect_size, num_loci, frac__mort_loci_dependent):
    num_ind = 6800
    max_age = 100
    mean_mort_rate = np.zeros((max_age))
    num_survivors = np.zeros((max_age))
    ages_of_interest = range(1,max_age)
    age_of_death_file = open('ages_of_death', 'w')
    for ind in range(num_ind):
        #age_freq_variation = 65 +  np.random.randint(2) * 20
        #p_locus_age = 1/ (age_freq_variation + 1) * np.ones(len(p_locus_age))
        simulation = sim_ind_binom(p_locus_age, locus_effect_size, ages_of_interest, num_loci, frac__mort_loci_dependent)
        age_at_death = simulation[0]
        ages_of_degredation = simulation[1]
        sim_mort_rates = simulation[2]
        for age in ages_of_interest:
            mean_mort_rate[age] += sim_mort_rates[age]
            if age_at_death > age:
                num_survivors[age]+=1
        print(age_at_death, file = age_of_death_file)
    filename = "age_structure" + str(num_causative) + "_" + str(int(len(p_locus_age) / sum(p_locus_age)))+ ".txt"

    age_struct = open(filename , "w")
    for age in ages_of_interest:
        if num_survivors[age] > 0:
            mean_mort_rate[age] = mean_mort_rate[age] / num_survivors[age]
            target_mort = math.exp(0.0807*age-9.575)
            print (age, num_survivors[age], mean_mort_rate[age], target_mort, sep = "\t", file = age_struct)
    return(mean_mort_rate)

def get_freqdeg_dist(param, num_loci, mean_p_locus_deg):
    freqdeg_dist = np.zeros((num_loci))
    mean_pdeg = 0
    for i in range(0, num_loci):
        #freqdeg_dist[i] = 10 * math.exp(param * (i / num_loci))
        freqdeg_dist[i] = param * (i / num_loci)
        #freqdeg_dist[i] = param
        mean_pdeg += 1 / freqdeg_dist[i]
    mean_pdeg = mean_pdeg / num_loci
    diff = (mean_pdeg - mean_p_locus_deg)**2
    print(param, mean_pdeg, mean_p_locus_deg, diff)
    return(diff)

def get_effect_size_dist(param, num_loci, num_causative, mean_effect_size):
    effect_size_dist = np.zeros((num_loci))
    mean_eff = 0
    for i in range(0,num_causative):
        #effect_size_dist[i] = mean_effect_size * max_locus_size_mult * math.exp(-1 * param * i / num_loci)
    #    if i / 2 == int(i/2):
        effect_size_dist[i] =  param
        mean_eff += effect_size_dist[i]
    diff = (mean_eff - mean_effect_size)**2
#    print(param, mean_eff, mean_effect_size* num_loci, diff)
    return(diff)


#num_causative_list = [301, 332, 371, 420, 484, 571, 697, 896, 1253, 2086, 6251]
num_causative_list = [301]
ages_of_interest = range(1,120)
mean_p_locus_deg =  0
p_locus_age = np.zeros((num_loci))
freqdeg_dist = np.zeros((num_loci))
#max_freqdeg = scipy.optimize.minimize_scalar(get_freqdeg_dist, bracket = (-20, 20), args = (num_loci, mean_p_locus_deg))
#print(max_freqdeg.x, mean_p_locus_deg)
#pdeg_param = scipy.optimize.minimize_scalar(get_freqdeg_dist, bracket = (0, 1000), args = (num_loci, mean_p_locus_deg))
test_variation_in_degredation_rate == 0
for replicate in range(0,1):
    for num_causative in num_causative_list:
        freq = 0.0807
        #freq = 0.0807 +  np.random.random() *0.04 - 0.02#75
        #freq = 75 +  np.random.random() *10 - 5#75
        for i in range(0, num_loci):
			if test_variation_in_degredation_rate ==1:
	            freqdeg_dist[i] = 1+199 * (i / num_loci)
	            p_locus_age[i] = 1/ (freqdeg_dist[i] + 1) #1 / freqdeg_dist[i]
	            p_locus_age[i] = 1/ (freq + 1) #1 / freqdeg_dist[i]
			else:
            	p_locus_age[i] = 1- math.exp(freq / (0.0807/math.log(1-1/76)))
            	mean_p_locus_deg += p_locus_age[i]
        mean_p_locus_deg = mean_p_locus_deg / num_loci
        print(mean_p_locus_deg, 1/mean_p_locus_deg)
        mean_effect_size = 1#-0.084 / (num_loci * math.log(1-mean_p_locus_deg) )

        #p_locus_age = np.ones((num_loci))*mean_p_locus_deg

        locus_effect_size = np.zeros((num_loci))

        #eff_size_param = scipy.optimize.minimize_scalar(get_effect_size_dist, bracket = (0, 1000), args = (num_loci, num_causative, mean_effect_size))
        #print(eff_size_param.x)
        for i in range(0, num_causative):
            locus_effect_size[i] = 1 / num_causative  #* i / 20000
        np.random.shuffle(locus_effect_size)
        np.random.shuffle(p_locus_age)

        frac__mort_loci_dependent = 0


        age_categories = range(10,90,10)
		#simulate_pop will generate a popualtion (6800 currently) and simulate then until they die. It exports the mean mortality rate and number of living individuals in 10 year bins
        simulate_pop(p_locus_age, locus_effect_size, num_loci, frac__mort_loci_dependent)
		#simulate cohorts simulates a cross sectional study; individuals are simulated and only included if they live until the target cohort age, as given by age_catagories.
		#simulate_cohorts exports the loci status of individuals at the time of sampling as well as each individuals age, mortality rate, and future age of death
        simulate_cohorts(p_locus_age, locus_effect_size, age_categories, frac__mort_loci_dependent, replicate)
