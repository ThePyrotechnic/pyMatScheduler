import copy
import random
from timeit import default_timer as timer
from collections import defaultdict, deque
from math import floor

real_names = [
    'Troy Lucas',
    'Charlene Olson',
    'Kathy Lewis',
    'Billie Campbell',
    'Emanuel Erickson',
    'Dale Singleton',
    'Orlando Shelton',
    'April Norris',
    'Craig Barber',
    'Shaun Ortiz',
    'Marilyn Sullivan',
    'Yvette Mcdaniel',
    'Penny Collier',
    'Oliver Walters',
    'Tyrone Lane',
    'Lawrence Johnston',
    'Brett Kelly',
    'Kathryn Morrison',
    'Brooke Thomas',
    'Bernard Haynes',
    'Nichole George',
    'Lorene Allison',
    'Agnes Mullins',
    'Roy Phelps',
    'Drew Obrien',
    'Gabriel Williams',
    'Opal Sherman',
    'Henrietta Frank',
    'Lauren Cortez',
    'Eddie Reese',
    'Hector Walsh',
    'Carroll Roberson',
    'Enrique Leonard',
    'Danny Hawkins',
    'Joey Robinson',
    'Sally Estrada',
    'Josefina Craig',
    'Hugo Mckenzie',
    'Eloise Gray',
    'Hazel Stevens',
    'Tom Harris',
    'Omar Franklin',
    'Luz Gilbert',
    'Philip Watts',
    'Eleanor Reed',
    'Dwight Hodges',
    'Alice Rogers',
    'Shirley Diaz',
    'Alexander Sutton',
    'Diana Guzman'
]

letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
           'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']


class Wrestler:
    """
    A generic wrestler, suitable for comparisons
    """
    age_tolerance = 3  # years
    weight_tolerance = 30  # pounds
    skill_tolerance = 3  # points

    def __init__(self, name: str, team: str, age: int = None, weight: int = None, skill: float = None):
        self.name = name
        self.team = team
        self.age = age if age is not None else random.randint(8, 15)
        self.skill = skill if skill is not None else random.uniform(1, 5)
        self.weight = weight if weight is not None else int(floor(random.randrange(95, 180, 5)))

    def __str__(self):
        # return '__\n{}\nTeam: {}\nAge: {}\nweight: {}\nskill: {}'.format(self.name, self.team, self.age, self.weight, round(self.skill, 1))
        return '{}'.format(self.name)

    def __repr__(self):
        # return "Wrestler(name = {}, team = {}, age = {}, weight = {}, skill = {})".format(self.name, self.team, self.age, self.weight, round(self.skill, 1))
        return '{}'.format(self.name)

    def __hash__(self):
        return hash((self.name, self.team, self.age, self.skill, self.weight))

    def __eq__(self, other):
        return (self.name, self.team, self.age, self.skill, self.weight) == (other.name, self.team, other.age, other.skill, other.weight)

    def __ne__(self, other):
        return not (self == other)

    def is_compatible(self, other: 'Wrestler') -> bool:
        return (other.age - self.age_tolerance < self.age < other.age + self.age_tolerance and
                other.weight - self.weight_tolerance < self.weight < other.weight + self.weight_tolerance and
                other.skill - self.skill_tolerance < self.skill < other.weight + self.skill_tolerance and
                self.team != other.team)


def generate_wrestlers(wrestler_count: int, names: list, teams: list) -> list:
    """
    Generates a list of random wrestlers.
    :param wrestler_count: The amount of wrestlers to generate
    :param names: A list of names to use
    :param teams: A list of teams to distribute evenly
    :return: A list of Wrestler objects
    """
    random.shuffle(names)
    return [Wrestler(name=name, team=teams[a % len(teams)]) for a, name in zip(range(wrestler_count), names)]


def calculate_matchups(wrestler_list: list) -> defaultdict:
    """
    Generates a dictionary with one entry for each wrestler.
    The value of each key is the list of wrestlers compatible with that wrestler.
    :param wrestler_list: A list of Wrestler objects
    :return: A dictionary, as described above
    """
    match_dict = defaultdict(list)
    for a, w in enumerate(wrestler_list):
        for cur in wrestler_list[a + 1:]:  # Don't need to start from the beginning each time.
            if w.is_compatible(cur):
                match_dict[w].append(cur)  # Add the possible match to each wrestler's dict (will not add repeats).
                match_dict[cur].append(w)

    return match_dict


def greatest_wrestler(match_dict: dict) -> Wrestler:
    """
    Finds the Wrestler with the most possible matchups in a given match dictionary.
    :param match_dict: The match dictionary to search
    :return: The Wrestler with the most possible matchups
    """
    return max(match_dict, key=lambda w: len(match_dict[w]))


def trim_matchups(match_dict: dict, limit: int):
    """
    Trims a match dictionary to fit the desired amount of rounds.
    Deletes matches between wrestlers with the most matches first, to ensure as many wrestlers get a match as possible.
    :param match_dict: The matchup dictionary to trim
    :param limit: The maximum amount of rounds to allow
    """
    w = greatest_wrestler(match_dict)
    while len(match_dict[w]) > limit:
        match_dict[w].sort(key=lambda cur_w: len(match_dict[cur_w]))
        greatest_match = match_dict[w].pop()
        match_dict[greatest_match].remove(w)
        w = greatest_wrestler(match_dict)


def remove_dupes(match_dict: dict):
    """
    Removes duplicate matches from the matchup dictionary
    :param match_dict: The dictionary to trim
    """
    for w, matchups in match_dict.items():
        for m in matchups:
            match_dict[m].remove(w)


def generate_match_queue(match_dict: dict, reverse: bool) -> deque:
    """
    Creates a queue of matches, with matchups of greatest combined match total at the front
    :param match_dict: The match dict to use
    :param reverse: Whether to reverse the order
    :return: A deque which satisfies the properties above
    """
    matches = []
    for w, matchups in match_dict.items():
        for m in matchups:
            matches.append((w, m, len(matchups) + len(match_dict[m])))

    q = deque(sorted(matches, key=lambda match: match[2], reverse=reverse))
    return q


def generate_schedule(match_dict: dict, mat_count: int, reverse: bool):
    """
    Generates a schedule of matches, given a matchup dictionary and a mat count
    :param match_dict: The matchup dictionary to use
    :param mat_count: The amount of mats (simultaneous matches)
    :param reverse: For debugging only, use False in production
    :return: A dictionary, where each key represents the schedule for a given mat
    """
    mat_dict = defaultdict(list)
    q = generate_match_queue(match_dict, reverse)
    cur_mat = 0
    first_match = None
    used_wrestlers = []
    while len(q) > 0:
        if cur_mat == mat_count or first_match == q[-1]:
            first_match = None
            cur_mat = 0
            q = generate_match_queue(match_dict, reverse)
            used_wrestlers = []
        match = q.pop()
        if not any(w in used_wrestlers for w in(match[0], match[1])):
            mat_dict[cur_mat].append((match[0], match[1]))
            cur_mat += 1
            used_wrestlers.extend((match[0], match[1]))
            match_dict[match[0]].remove(match[1])
        else:
            q.appendleft(match)
            if first_match is None:
                first_match = match

    return mat_dict


def main(num_wrestlers: int, num_teams: int, mat_count: int, max_matches_per_wrestler: int):
    avg_n_spread, avg_rev_spread, avg_n_max, avg_rev_max, avg_n_min, avg_rev_min = 0, 0, 0, 0, 0, 0
    n = 100

    numbers = ['{}'.format(n) for n in range(num_wrestlers)]
    teams = ['{}'.format(n) for n in range(num_teams)]

    start = timer()
    for a in range(n):
        match_dict = calculate_matchups(generate_wrestlers(num_wrestlers, numbers, teams))
        # match_dict = {
        #     'T': ['E'],
        #     'E': ['T', 'B', 'G', 'M', 'J'],
        #     'V': ['G', 'M'],
        #     'G': ['V', 'E', 'U'],
        #     'M': ['V', 'E', 'U'],
        #     'B': ['E'],
        #     'U': ['G', 'M'],
        #     'J': ['E']
        # }

        trim_matchups(match_dict, max_matches_per_wrestler)
        remove_dupes(match_dict)
        rev_dict = copy.deepcopy(match_dict)
        raw_schedule = generate_schedule(match_dict, mat_count, False)
        rev_raw_schedule = generate_schedule(rev_dict, mat_count, True)

        cur_max = len(raw_schedule[0])
        cur_min = len(raw_schedule[2])
        avg_n_spread += cur_max - cur_min
        avg_n_max += cur_max
        avg_n_min += cur_min

        cur_max = len(rev_raw_schedule[0])
        cur_min = len(rev_raw_schedule[2])
        avg_rev_spread += cur_max - cur_min
        avg_rev_max += cur_max
        avg_rev_min += cur_min

    end = timer()
    avg_n_spread /= n
    avg_n_max /= n
    avg_n_min /= n
    avg_rev_spread /= n
    avg_rev_max /= n
    avg_rev_min /= n

    print('{} test averages:\nn_spread: \t\t{}|n_max: \t{}|n_min \t{}\nrev_spread: \t{}|rev_max: \t{}|rev_min \t{}'.format(
        n, avg_n_spread, avg_n_max, avg_n_min, avg_rev_spread, avg_rev_max, avg_rev_min))
    print('Total time: {}'.format(end - start))


if __name__ == '__main__':
    # num_wrestlers, num_teams, mat_count, max_matches_per_wrestler
    main(210, 3, 3, 4)

    # T: E
    # E: T B G M J
    # V: G M
    # G: V E U
    # M: V E U
    # B: E
    # U: G M
    # J: E

    # Good test case to check by hand
    # match_dict = {
    #     'F': ['C', 'D', 'A', 'E'],
    #     'A': ['C', 'D', 'F'],
    #     'E': ['F'],
    #     'C': ['A', 'F'],
    #     'B': ['G'],
    #
    #     'G': ['B', 'D'],
    #     'D': ['A', 'F', 'G']
    # }
