from graph import Network
import numpy as np
import logging
from sys import stdout
import gc

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import settings

formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
handler = logging.StreamHandler(stdout)
handler.setFormatter(formatter)
logger = logging.getLogger('AGDN')
logger.setLevel(logging.INFO)
logger.addHandler(handler)


def main():
    population = 9e6
    ventilators = 2e3
    days = np.linspace(0, 364, 365, endpoint=True)
    i_peak = np.array([ventilators for _ in days])
    for eta in settings.ETAS:
        logger.info("initializing network with eta {}".format(eta))
        network = Network.initialize(zones=260, total_population=population, eta=eta)
        infectious = list()
        for day in days:
            logger.info("starting day {}".format(day))
            infectious.append(sum(map(lambda zone: zone.population.hospitalized, network.vertices)))
            network.increment()
        logger.info("objects in mem {}".format(gc.get_count()))
        gc.collect()
        plt.plot(days, infectious)

    plt.plot(days, i_peak)
    legends = ['H(t); eta = {}'.format(eta) for eta in settings.ETAS]
    legends.append('I peak')
    plt.legend(legends)
    plt.show()


if __name__ == '__main__':
    main()
