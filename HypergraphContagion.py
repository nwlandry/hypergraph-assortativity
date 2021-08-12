import Hypergraph
import numpy as np
import random
import multiprocessing as mp
from collections import defaultdict
from collections import Counter


def getEquilibriumPoint(G, gamma, beta, fraction_infected, tmax, fractionToAverage, numSimulations, isVerbose):
    equilibriumPoint = 0
    for sim in range(numSimulations):
        t, S, I = Gillespie_SIS(G, beta, gamma, transmission_function=individual_contagion, rho=fraction_infected, tmin=0, tmax=tmax)
        equilibriumPoint += averageEquilibrium(t, I, averagingTime=fractionToAverage*(np.max(t)-np.min(t)))/numSimulations
    if isVerbose:
        print(equilibriumPoint, flush=True)
    return equilibriumPoint

def averageEquilibrium(time, data, averagingTime=None):
    n = np.size(time, axis=0)
    if averagingTime is not None:
        cutoff = time[-1] - averagingTime
        i = n - 1
        while cutoff < time[i] and i > 0:
            i -= 1
        # interpolate
        if i > -np.size(time, axis=0):
            firstDataValue = np.interp(cutoff, [time[i], time[i+1]], [data[i], data[i+1]])
            time = np.concatenate([np.array([cutoff]), time[i+1:]])
            data = np.concatenate([np.array([firstDataValue]), data[i+1:]])
    try:
        return np.average(data, weights=[(np.max(time)-np.min(time))/np.size(time, axis=0)] + list(time[1:] - time[:-1]))
    except:
        return np.mean(data)

#######################
#                     #
#   Auxiliary stuff   #
#                     #
#######################

# built-in functions
def collective_contagion(status, neighbors):
    for i in neighbors:
        if status[i] != 'I':
            return 0
    return 1

def individual_contagion(status, neighbors):
    for i in neighbors:
        if status[i] == 'I':
            return 1
    return 0

def threshold(status, neighbors, threshold=0.5):
    meanContagion = sum([status[i] == 'I' for i in neighbors])/len(neighbors)
    if meanContagion < threshold:
        return 0
    elif meanContagion >= threshold:
        return 1

def majority_vote(status, neighbors):
    meanContagion = sum([status[i] == 'I' for i in neighbors])/len(neighbors)
    if meanContagion < 0.5:
        return 0
    elif meanContagion > 0.5:
        return 1
    else:
        return random.choice([0, 1])

def size_dependent(status, neighbors):
    return sum([status[i] == 'I' for i in neighbors])

# distribution function

def _truncated_exponential_(rate, T):
    r'''returns a number between 0 and T from an
    exponential distribution conditional on the outcome being between 0 and T'''
    t = random.expovariate(rate)
    L = int(t/T)
    return t - L*T

# Classes for the simulations

class _ListDict_(object):
    def __init__(self, weighted = False):
        self.item_to_position = {}
        self.items = []

        self.weighted = weighted
        if self.weighted:
            self.weight = defaultdict(int) #presume all weights positive
            self.max_weight = 0
            self._total_weight = 0
            self.max_weight_count = 0


    def __len__(self):
        return len(self.items)

    def __contains__(self, item):
        return item in self.item_to_position

    def _update_max_weight(self):
        C = Counter(self.weight.values())  #may be a faster way to do this, we only need to count the max.
        self.max_weight = max(C.keys())
        self.max_weight_count = C[self.max_weight]


    def insert(self, item, weight = None):
        r'''
        If not present, then inserts the thing (with weight if appropriate)
        if already there, replaces the weight unless weight is 0

        If weight is 0, then it removes the item and doesn't replace.

        WARNING:
            replaces weight if already present, does not increment weight.


        '''
        if self.__contains__(item):
            self.remove(item)
        if weight != 0:
            self.update(item, weight_increment=weight)


    def update(self, item, weight_increment = None):
        r'''
        If not present, then inserts the thing (with weight if appropriate)
        if already there, increments weight

        WARNING:
            increments weight if already present, cannot overwrite weight.
        '''
        if weight_increment is not None: #will break if passing a weight to unweighted case
            if weight_increment > 0 or self.weight[item] != self.max_weight:
                self.weight[item] = self.weight[item] + weight_increment
                self._total_weight += weight_increment
                if self.weight[item] > self.max_weight:
                    self.max_weight_count = 1
                    self.max_weight = self.weight[item]
                elif self.weight[item] == self.max_weight:
                    self.max_weight_count += 1
            else: #it's a negative increment and was at max
                self.max_weight_count -= 1
                self.weight[item] = self.weight[item] + weight_increment
                self._total_weight += weight_increment
                self.max_weight_count -= 1
                if self.max_weight_count == 0:
                    self._update_max_weight
        elif self.weighted:
            raise Exception('if weighted, must assign weight_increment')

        if item in self: #we've already got it, do nothing else
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1

    def remove(self, choice):
        position = self.item_to_position.pop(choice) # why don't we pop off the last item and put it in the choice index?
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position

        if self.weighted:
            weight = self.weight.pop(choice)
            self._total_weight -= weight
            if weight == self.max_weight:
                self.max_weight_count -= 1
                if self.max_weight_count == 0 and len(self)>0:
                    self._update_max_weight()

    def choose_random(self):
        if self.weighted:
            while True:
                choice = random.choice(self.items)
                if random.random() < self.weight[choice]/self.max_weight:
                    break
            return choice

        else:
            return random.choice(self.items)


    def random_removal(self):
        r'''uses other class methods to choose and then remove a random node'''
        choice = self.choose_random()
        self.remove(choice)
        return choice

    def total_weight(self):
        if self.weighted:
            return self._total_weight
        else:
            return len(self)
    def update_total_weight(self):
        self._total_weight = sum(self.weight[item] for item in self.items)

##########################
#                        #
#    SIMULATION CODE     #
#                        #
##########################

'''
    The code in the region below is used for stochastic simulation of
    epidemics on networks
'''
def Gillespie_SIS_perturbed(G, tau, gamma, transmission_function=collective_contagion, initial_infecteds=None, rho=None, tmin=0, tmax=100, recovery_weight=None, transmission_weight=None, **args):

    if transmission_weight is not None:
        def edgeweight(item):
            return item[transmission_weight]
    else:
        def edgeweight(item):
            return None

    if recovery_weight is not None:
        def nodeweight(u):
            return G.nodes[u][recovery_weight]
    else:
        def nodeweight(u):
            return None

    gamma = float(gamma)

    if initial_infecteds is None:
        if rho is None:
            initial_number = 1
        else:
            initial_number = int(round(G.number_of_nodes()*rho))
        initial_infecteds=random.sample(G.nodeLabels, initial_number)
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds] # for a single node

    I = [len(initial_infecteds)]
    S = [G.number_of_nodes()-I[0]]
    times = [tmin]

    t = tmin
    transmissions = []

    status = defaultdict(lambda : 'S')
    for node in initial_infecteds:
        status[node] = 'I'

    if recovery_weight is None:
        infecteds = _ListDict_()
    else:
        infecteds = _ListDict_(weighted=True)

    IS_links = dict()
    for size in G.getHyperedgeSizes():
        if transmission_weight is None:
            IS_links[size] = _ListDict_()
        else:
            IS_links[size] = _ListDict_(weighted=True)

    for node in initial_infecteds:
        infecteds.update(node, weight_increment = nodeweight(node))
        for uid, nbrData in G.neighbors[node].items():  #must have this in a separate loop after assigning status of node
            for nbr in nbrData["neighbors"]: # there may be self-loops so account for this later
                if status[nbr] == 'S':
                    contagion = transmission_function(status, G.neighbors[nbr][uid]["neighbors"], **args)
                    if contagion != 0:
                        IS_links[len(nbrData["neighbors"])+1].update((G.neighbors[nbr][uid]["neighbors"], nbr), weight_increment=edgeweight(nbrData)) # need to be able to multiply by the contagion?

    total_rates = dict()
    total_rates[1] = gamma*infecteds.total_weight()#I_weight_sum
    for size in G.getHyperedgeSizes():
        total_rates[size] = tau[size]*IS_links[size].total_weight() #IS_weight_sum

    total_rate = sum(total_rates.values())

    delay = random.expovariate(total_rate)
    t += delay

    while t <= tmax:
        # rejection sampling
        while True:
            choice = random.choice(list(total_rates.keys()))
            if random.random() < total_rates[choice]/total_rate:
                break

        if choice == 1: #recover
            if I[-1] > 1:
                recovering_node = infecteds.random_removal() # chooses a node at random and removes it
                status[recovering_node] = 'S' # update node status

                # Find the SI links for the recovered node to get reinfected
                for uid, nbrData in G.neighbors[recovering_node].items():
                    contagion =  transmission_function(status, nbrData["neighbors"], **args)
                    if contagion != 0:
                        IS_links[len(nbrData["neighbors"])+1].update((nbrData["neighbors"], recovering_node), weight_increment = edgeweight(nbrData))

                # reduce the number of infected links because of the healing
                for uid, nbrData in G.neighbors[recovering_node].items():
                    for nbr in nbrData["neighbors"]:
                        if status[nbr] == 'S' and (G.neighbors[nbr][uid]["neighbors"], nbr) in IS_links[len(nbrData["neighbors"])+1]: # if the key doesn't exist, don't attempt to remove it
                            contagion = transmission_function(status, G.neighbors[nbr][uid]["neighbors"], **args)
                            if contagion == 0:
                                try:
                                    IS_links[len(nbrData["neighbors"])+1].remove((G.neighbors[nbr][uid]["neighbors"], nbr))
                                except:
                                    pass
                times.append(t)
                S.append(S[-1]+1)
                I.append(I[-1]-1)
            else:
                times.append(t)
                S.append(S[-1])
                I.append(I[-1])
        else:
            transmitter, recipient = IS_links[choice].choose_random()
            status[recipient]='I'

            infecteds.update(recipient, weight_increment = nodeweight(recipient))

            for uid, nbrData in G.neighbors[recipient].items():
                try:
                    IS_links[len(nbrData["neighbors"])+1].remove((nbrData["neighbors"], recipient)) # multiply by contagion?
                except:
                    pass

            for uid, nbrData in G.neighbors[recipient].items():
                for nbr in nbrData["neighbors"]:
                    if status[nbr] == 'S':
                        contagion = transmission_function(status, G.neighbors[nbr][uid]["neighbors"], **args)
                        if contagion != 0:
                            IS_links[len(nbrData["neighbors"])+1].update((G.neighbors[nbr][uid]["neighbors"], nbr), weight_increment = edgeweight(nbrData))

            times.append(t)
            S.append(S[-1]-1)
            I.append(I[-1]+1)

        total_rates[1] = gamma*infecteds.total_weight()#I_weight_sum
        for size in G.getHyperedgeSizes():
            total_rates[size] = tau[size]*IS_links[size].total_weight() #IS_weight_sum
        total_rate = sum(total_rates.values())
        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
    return np.array(times), np.array(S), np.array(I)


def Gillespie_SIS(G, tau, gamma, transmission_function=collective_contagion, initial_infecteds=None, rho=None, tmin=0, tmax=100, recovery_weight=None, transmission_weight=None, **args):

    if transmission_weight is not None:
        def edgeweight(item):
            return item[transmission_weight]
    else:
        def edgeweight(item):
            return None

    if recovery_weight is not None:
        def nodeweight(u):
            return G.nodes[u][recovery_weight]
    else:
        def nodeweight(u):
            return None

    gamma = float(gamma)

    if initial_infecteds is None:
        if rho is None:
            initial_number = 1
        else:
            initial_number = int(round(G.number_of_nodes()*rho))
        initial_infecteds=random.sample(G.nodeLabels, initial_number)
    elif G.has_node(initial_infecteds):
        initial_infecteds=[initial_infecteds] # for a single node

    I = [len(initial_infecteds)]
    S = [G.number_of_nodes()-I[0]]
    times = [tmin]

    t = tmin
    transmissions = []

    status = defaultdict(lambda : 'S')
    for node in initial_infecteds:
        status[node] = 'I'

    if recovery_weight is None:
        infecteds = _ListDict_()
    else:
        infecteds = _ListDict_(weighted=True)

    IS_links = dict()
    for size in G.getHyperedgeSizes():
        if transmission_weight is None:
            IS_links[size] = _ListDict_()
        else:
            IS_links[size] = _ListDict_(weighted=True)

    for node in initial_infecteds:
        infecteds.update(node, weight_increment = nodeweight(node))
        for uid, nbrData in G.neighbors[node].items():  #must have this in a separate loop after assigning status of node
            for nbr in nbrData["neighbors"]: # there may be self-loops so account for this later
                if status[nbr] == 'S':
                    contagion = transmission_function(status, G.neighbors[nbr][uid]["neighbors"], **args)
                    if contagion != 0:
                        IS_links[len(nbrData["neighbors"])+1].update((G.neighbors[nbr][uid]["neighbors"], nbr), weight_increment=edgeweight(nbrData)) # need to be able to multiply by the contagion?

    total_rates = dict()
    total_rates[1] = gamma*infecteds.total_weight()#I_weight_sum
    for size in G.getHyperedgeSizes():
        total_rates[size] = tau[size]*IS_links[size].total_weight() #IS_weight_sum

    total_rate = sum(total_rates.values())

    delay = random.expovariate(total_rate)
    t += delay

    while infecteds and t < tmax:
        # rejection sampling
        while True:
            choice = random.choice(list(total_rates.keys()))
            if random.random() < total_rates[choice]/total_rate:
                break
        #choice = np.random.choice(total_rates.keys(), p=np.array(total_rates.values())/sum(total_rates.values()))
        #choice = random.choices(list(total_rates.keys()), weights=list(total_rates.values()), k=1)[0]

        if choice == 1: #recover
            recovering_node = infecteds.random_removal() # chooses a node at random and removes it
            status[recovering_node] = 'S' # update node status

            # Find the SI links for the recovered node to get reinfected
            for uid, nbrData in G.neighbors[recovering_node].items():
                contagion =  transmission_function(status, nbrData["neighbors"], **args)
                if contagion != 0:
                    IS_links[len(nbrData["neighbors"])+1].update((nbrData["neighbors"], recovering_node), weight_increment = edgeweight(nbrData))

            # reduce the number of infected links because of the healing
            for uid, nbrData in G.neighbors[recovering_node].items():
                for nbr in nbrData["neighbors"]:
                    if status[nbr] == 'S' and (G.neighbors[nbr][uid]["neighbors"], nbr) in IS_links[len(nbrData["neighbors"])+1]: # if the key doesn't exist, don't attempt to remove it
                        contagion = transmission_function(status, G.neighbors[nbr][uid]["neighbors"], **args)
                        if contagion == 0:
                            try:
                                IS_links[len(nbrData["neighbors"])+1].remove((G.neighbors[nbr][uid]["neighbors"], nbr)) # should this be "update" instead of "remove"?
                            except:
                                pass

            times.append(t)
            S.append(S[-1]+1)
            I.append(I[-1]-1)
        else:
            transmitter, recipient = IS_links[choice].choose_random()
            status[recipient]='I'

            infecteds.update(recipient, weight_increment = nodeweight(recipient))

            for uid, nbrData in G.neighbors[recipient].items():
                try:
                    IS_links[len(nbrData["neighbors"])+1].remove((nbrData["neighbors"], recipient)) # multiply by contagion?
                except:
                    pass

            for uid, nbrData in G.neighbors[recipient].items():
                for nbr in nbrData["neighbors"]:
                    if status[nbr] == 'S':
                        contagion = transmission_function(status, G.neighbors[nbr][uid]["neighbors"], **args)
                        if contagion != 0:
                            IS_links[len(nbrData["neighbors"])+1].update((G.neighbors[nbr][uid]["neighbors"], nbr), weight_increment = edgeweight(nbrData))
            times.append(t)
            S.append(S[-1]-1)
            I.append(I[-1]+1)

        total_rates[1] = gamma*infecteds.total_weight()#I_weight_sum
        for size in G.getHyperedgeSizes():
            total_rates[size] = tau[size]*IS_links[size].total_weight() #IS_weight_sum
        total_rate = sum(total_rates.values())
        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
    return np.array(times), np.array(S), np.array(I)