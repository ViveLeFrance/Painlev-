from dataclasses import dataclass
from typing import List, Set, Tuple



class Crossing:
    def __init__(self, name: str):
        self.name = name

    def __hash__(self):
        return hash(self.name)
    
    def __repr__(self):
        return self.name
    

PLUS = '+'
MINUS = '-'


class Sector:
    def __init__(self, sign: str):
        self.sign = sign
        
    def __hash__(self):
        return hash(self.sign)
    
    def __repr__(self) -> str:
        return self.sign

plus = Sector(PLUS)
minus = Sector(MINUS)


class Loop:
    def __init__(self, crossing: Crossing):
        self.crossing = crossing

    def __hash__(self):
        return hash(self.crossing)
    
    def __repr__(self):
        return f'loop_{self.crossing.name}'


class Move:
    def __init__(self, crossing: Crossing, sector: Sector, loops: List[Loop]):
        self.crossing = crossing
        self.sector = sector
        self.loops = loops

    def __hash__(self):
        return hash((self.crossing, self.sector, tuple(self.loops)))

    def __repr__(self):
        return(f'Move({self.crossing.__repr__()}, {self.sector.__repr__()}, {self.loops.__repr__()}')



x = Crossing('x')
y = Crossing('y')
z = Crossing('z')
p = Crossing('p')
q = Crossing('q')
r = Crossing('r')
s = Crossing('s')
a = Crossing('a')
b = Crossing('b')
c = Crossing('c')

l_a = Loop(a)
l_b = Loop(b)
l_c = Loop(c)


basepoints = {
    't_1': (Move())
}

legal_moves: dict[Tuple, Move] = {
    (x,plus): [
        Move(y, plus, []),
        Move(a, minus, [l_a]),
        Move(x, plus, [l_a]),
        Move(r, plus, [l_a]),
        Move(z, plus, [l_a]),
        Move(b, plus, [l_a, l_b]),
        ],
    (x,minus): [
        Move(z, minus, [l_b]),
        Move(r, minus, [l_b])
        ],
    (y,plus): [
        Move(r, plus, []),
        Move(z, plus, []),
        Move(b, plus, [])
    ],
    (y, minus): [
        Move(x, minus, [])
        ],
    (z, plus): [
        Move(a, minus, [l_a]),
        Move(x, plus, [l_a]),
        Move(r, plus, [l_a]),
        Move(b, plus, [l_a, l_b]),
        ],
    (z, minus): [
        Move(y, minus, []),
        Move(z, minus, [l_b]),
        Move(b, minus, [l_b]),
        Move(r, minus, [l_b]),
        Move(x, minus, [l_b]),
    ],
    (p, plus): [],
    (p, minus): [
        Move(x, plus, []),
        Move(r, plus, []),
        Move(z, plus, []),
        Move(b, plus, [l_b])
    ],
    (q, plus): [
        Move(r, plus, []),
        Move(b, plus, [l_b])
    ],
    (q, minus): [
        Move(x, minus, []),
    ],
    (r, plus): [],
    (r, minus): [
        Move(q, minus, []),
        Move(p, minus, []),
        Move(c, plus, [])
    ],
    (s, plus): [],
    (s, minus): [
        Move(q, minus, []),
        Move(p, minus, []),
        Move(c, plus, [])
    ],
    (a, plus): [
        Move(a, minus, [l_a]),
        Move(s, minus, [l_a]),
        Move(x, plus, [l_a]),
        Move(r, plus, [l_a]),
        Move(z, plus, [l_a]),
        Move(b, plus, [l_a, l_b])
    ],
    (a, minus): [
        Move(y, minus, []),
        Move(b, minus, [l_b]),
        Move(z, minus, [l_b]),
        Move(r, minus, [l_b]),
        Move(x, minus, [l_b])
    ],
    (b, plus): [
        Move(b, minus, [l_b]),
        Move(z, minus, [l_b]),
        Move(r, minus, [l_b]),
        Move(x, minus, [l_b])
    ],
    (b, minus): [
        Move(y, plus, []),
        Move(a, minus, [l_a]),
        Move(x, plus, [l_a]),
        Move(r, plus, [l_a]),
        Move(z, plus, [l_a]),
        Move(b, plus, [l_a, l_b])
    ],
    (c, plus): [],
    (c, minus): []
    
}

start_moves = [
    Move(q, minus, []),
    Move(p, minus, []),
    Move(c, plus, [])
]


def generate_permissible_words(legal_moves, start_moves, max_depth):
    """
    Generates all permissible words based on the start moves, legal moves, and loop constraints.

    Parameters:
    - legal_moves: The dictionary of legal moves.
    - start_moves: A list of Move objects where the word can start.
    - max_depth: The maximum length of the generated word.

    Returns:
    - A list of permissible sequences of moves that start from a given set of start moves and end at crossing 'c'.
    """

    def traverse(current_crossing, current_sector, current_word, traversed_loops: Set[Loop], depth):
        # Base case: If the maximum depth is reached, return the word only if it ends at crossing 'c'
        if depth == max_depth or current_crossing == c:
            return [current_word]
        

        words = []
        # Get the legal moves from the current crossing-sector pair
        for move in legal_moves.get((current_crossing, current_sector)):
            assert isinstance(move, Move)
            # Check if any loop in the current move has already been traversed
            if all(loop not in traversed_loops for loop in move.loops):
                # Update the current word and traversed loops
                new_word = current_word + [move.crossing]
                new_traversed_loops = traversed_loops | set(move.loops)
                # Recursively traverse to find all permissible words
                words.extend(traverse(move.crossing, move.sector, new_word, new_traversed_loops, depth + 1))

        return words

    # Initialize the traversal from all start moves
    all_words = []
    for move in start_moves:
        # Start traversing from each of the start moves
        words_from_start = traverse(move.crossing, move.sector, [move.crossing], set(), 1)
        all_words.extend(words_from_start)

    return all_words



perm_words = generate_permissible_words(legal_moves, start_moves, 50)
for word in perm_words:
    print(word)