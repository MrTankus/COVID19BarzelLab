import itertools
import numpy as np
from scipy.integrate import odeint

import settings


class State(object):
    RED = 0
    YELLOW = 1
    GREEN = 2


class ZonePopulation(object):

    def __init__(self, susceptible, exposed, infected, recovered, hospitalized):
        self.susceptible = susceptible
        self.exposed = exposed
        self.infected = infected
        self.hospitalized = hospitalized
        self.recovered = recovered
        self.total = susceptible + exposed + infected + recovered

    def calculate_positive_fraction(self):
        return (self.exposed + self.infected + self.hospitalized) / self.total


class MigrationPopulation(object):

    def __init__(self, fraction):
        assert 0 < fraction < 1
        self.fraction = fraction

    @classmethod
    def calculate_migration_coefficient(cls, from_zone, to_zone, all_zones):
        r = from_zone.distance_from(to_zone)
        rhs = from_zone.population.total * sum([zone.population.total * r ** -2 for zone in filter(lambda z: z != from_zone, all_zones)])
        return from_zone.population.total * settings.MIGRATION_COEFICIENT_PERCENTAGE / rhs

    def migrate_population(self, from_population, to_population):

        def normalize(n):
            if n > 0:
                return int(n * self.fraction) if n > 10 else 1
            else:
                return 0
        migrated_susceptible = normalize(from_population.susceptible)
        migrated_exposed = normalize(from_population.exposed)
        migrated_infected = 0
        migrated_hospitalized = 0
        migrated_recovered = normalize(from_population.recovered)

        to_population.susceptible += migrated_susceptible
        to_population.exposed += migrated_exposed
        to_population.infected += migrated_infected
        to_population.hospitalized += migrated_hospitalized
        to_population.recovered += migrated_recovered

        from_population.susceptible -= migrated_susceptible
        from_population.exposed -= migrated_exposed
        from_population.infected -= migrated_infected
        from_population.hospitalized -= migrated_hospitalized
        from_population.recovered -= migrated_recovered


class Zone(object):

    def __init__(self, id, location, population, state=None):
        self.id = id
        self.location = location
        self.population = population
        self.state = state
        self.migration_to_map = dict()
        self.migration_from_map = dict()

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return self.id == other.id

    def distance_from(self, vertex):
        return np.sqrt((self.location[0] - vertex.location[0]) ** 2 + (self.location[1] - vertex.location[1]) ** 2)

    def run_internal_dynamics(self, alpha, beta, gamma, delta):

        def model(z, t, a, b, g, d, population):
            s = z[0]
            e = z[1]
            i = z[2]
            h = z[3]
            dsdt = -b * s * (e / population)
            dedt = b * e * (s / population) - g * e
            didt = g * e - a * i
            dhdt = d * i - a * h
            drdt = a * (i + h)
            return [dsdt, dedt, didt, dhdt, drdt]

        t = np.linspace(0, 9, 10 * 60)
        state = [self.population.susceptible, self.population.exposed, self.population.infected, self.population.hospitalized, self.population.recovered]
        z = odeint(model, state, t, args=(alpha, beta, gamma, delta, self.population.total))
        susceptible_over_time = z[:, 0]
        exposed_over_time = z[:, 1]
        infected_over_time = z[:, 2]
        hospitalzied_over_time = z[:, 2]
        recovered_over_time = z[:, 3]

        self.population.susceptible = susceptible_over_time[-1]
        self.population.exposed = exposed_over_time[-1]
        self.population.infected = infected_over_time[-1]
        self.population.hospitalized = hospitalzied_over_time[-1]
        self.population.recovered = recovered_over_time[-1]

    def set_state(self, eta):
        calc = (self.population.exposed + self.population.infected) / self.population.total
        if calc < eta:
            self.state = State.GREEN
        else:
            self.state = State.RED


class Edge(object):

    def __init__(self, v1, v2, weight):
        self.v1 = v1
        self.v2 = v2
        self.weight = weight

    def is_on_edge(self, vertex):
        return self.v1 == vertex or self.v2 == vertex


class Network(object):

    def __init__(self, vertices, edges, eta):
        self.vertices = set(vertices)
        self.edges = set(edges)
        self.vertices_map = dict([(v.id, v) for v in self.vertices])
        self.eta = eta

    @classmethod
    def initialize(cls, zones, total_population, eta):
        vertices = set()
        i = 0
        j = 0
        for zone in range(zones):
            # TODO - get better way to randomize zones populations
            zone_population = int(total_population / zones)
            if zone_population == 0:
                zone_population = total_population
            else:
                total_population -= zone_population

            # TODO - randomize local initial state.

            population = ZonePopulation(susceptible=zone_population - 1, exposed=1,
                                        infected=0, hospitalized=0, recovered=0)
            zone = Zone(id=zone, location=(i, j), population=population, state=State.GREEN)
            j += 1
            if j % 3 == 0:
                i += 1
                j = 0
            vertices.add(zone)
        edges = itertools.combinations(vertices, 2)

        edges = set([Edge(v1=p[0], v2=p[1], weight=((p[0].population.total * p[1].population.total) * (p[0].distance_from(p[1]) ** -2))) for p in edges])

        for from_zone in vertices:
            from_zone_edges = filter(lambda e: e.is_on_edge(from_zone), edges)
            for link in from_zone_edges:
                to_zone = link.v1 if from_zone != link.v1 else link.v2
                c = MigrationPopulation.calculate_migration_coefficient(from_zone=from_zone, to_zone=to_zone, all_zones=vertices)
                # TODO - fix this! link.weight is now c * (p1 * p2 / r^2) and not (1/r^2)
                if np.random.random() < c * link.weight:
                    migration_population = MigrationPopulation(fraction=settings.MIGRATION_PERCENTAGE)
                    from_zone.migration_to_map[to_zone.id] = migration_population
                    to_zone.migration_from_map[from_zone.id] = migration_population

        return Network(vertices=vertices, edges=edges, eta=eta)

    def increment(self):
        for zone in self.vertices:
            zone.set_state(eta=self.eta)
        for zone in self.vertices:
            self.migrate_to(zone)
        for zone in self.vertices:
            if zone.state is not None:
                if zone.state != State.RED:
                    zone.run_internal_dynamics(alpha=settings.ALPHA, beta=settings.BETA, gamma=settings.GAMMA, delta=settings.DELTA)
                else:
                    zone.run_internal_dynamics(alpha=settings.ALPHA, beta=0, gamma=settings.GAMMA, delta=settings.DELTA)
        for zone in self.vertices:
            self.migrate_back(zone)

    def migrate_to(self, from_zone):
        for to_zone_id in from_zone.migration_to_map:
            to_zone = self.vertices_map.get(to_zone_id)
            migration_population = from_zone.migration_to_map.get(to_zone_id)
            if to_zone:
                migration_population.migrate_population(from_population=from_zone.population, to_population=to_zone.population)

    def migrate_back(self, to_zone):
        for from_zone_id in to_zone.migration_from_map:
            from_zone = self.vertices_map.get(from_zone_id)
            migration_population = to_zone.migration_from_map.get(from_zone_id)
            if from_zone:
                migration_population.migrate_population(from_population=from_zone.population,
                                                        to_population=to_zone.population)
