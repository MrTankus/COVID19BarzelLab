from graph import Network
import numpy as np
import logging
from sys import stdout
import gc

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
handler = logging.StreamHandler(stdout)
handler.setFormatter(formatter)
logger = logging.getLogger('AGDN')
logger.setLevel(logging.INFO)
logger.addHandler(handler)


def main():
    population = 9e6
    ventilators = 2e3
    k = 0.1
    days = np.linspace(0, 364, 365, endpoint=True)
    i_peak = np.array([ventilators / k for _ in days])
    etas = [0.01, 0.001]
    for eta in etas:
        logger.info("initializing network with eta {}".format(eta))
        network = Network.initialize(zones=260, total_population=population, eta=eta)
        infectious = list()
        for day in days:
            logger.info("starting day {}".format(day))
            infectious.append(sum(map(lambda zone: zone.population.infected, network.vertices)))
            network.increment()
        logger.info("objects in mem {}".format(gc.get_count()))
        gc.collect()
        plt.plot(days, infectious)

    plt.plot(days, i_peak)
    plt.legend(['infectious over days - 0.01', 'infectious over days - 0.001', 'I peak'])
    plt.show()


if __name__ == '__main__':
    main()
