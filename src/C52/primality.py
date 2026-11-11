"""implements probabilistic primality tests"""

from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from timeit import repeat
from random import Random
from collections import deque
from math import gcd, isqrt


OPTIMAL_BENCHMARK_REPEAT_AMT = 5
DEFAULT_BENCHMARK_ITER_AMT = 1000
DEFAULT_ROUND_AMT = 5

MAX_WITNESSES = 5

MS_IN_S = 1000


class PrimalityTestABC(ABC):
    rounds: int
    
    def __init__(self, rounds: int = DEFAULT_ROUND_AMT) -> None:
        if rounds < 1:
            raise ValueError('Number of rounds must be greater than 0')

        self.rounds = rounds

    @property
    @abstractmethod
    def algorithm_name(self) -> str:
        pass

    def benchmark_s(self, n: int, iter: int = DEFAULT_BENCHMARK_ITER_AMT) -> float:
        return min(repeat(
            lambda: self.is_prime(n),
            repeat=OPTIMAL_BENCHMARK_REPEAT_AMT,
            number=iter
        ))

    def benchmark_ms(self, n: int, iter: int = DEFAULT_BENCHMARK_ITER_AMT) -> float:
        return self.benchmark_s(n, iter) * MS_IN_S

    def gather_stats(self, n: int, iter: int = DEFAULT_BENCHMARK_ITER_AMT) -> int:
        res = 0

        for _ in range(iter):
            if self.is_prime(n).is_prime:
                res += 1

        return res

    def is_prime(self, n: int) -> TestResult:
        if n < 2:
            return TestResult(self.algorithm_name, n, is_prime=False, reason='Numbers less than 2 are not prime')

        if n in (2, 3):
            return TestResult(self.algorithm_name, n, reason='2 and 3 are prime')

        if n % 2 == 0:
            return TestResult(self.algorithm_name, n, is_prime=False, reason='Even numbers greater than 2 are composite')

        return self._run_test(n)

    def _run_test(self, n: int) -> TestResult:
        rng = Random()
        res = TestResult(self.algorithm_name, n)

        for _ in range(self.rounds):
            res.rounds_executed += 1

            if not self._perform_round(rng, res):
                break

        return res

    @abstractmethod
    def _perform_round(self, rng: Random, res: TestResult) -> bool:
        pass


class FermatTest(PrimalityTestABC):
    @property
    def algorithm_name(self) -> str:
        return 'Fermat'
    
    def _perform_round(self, rng: Random, res: TestResult) -> bool:
        n = res.candidate
        a = rng.randint(2, n - 2)

        res.witnesses.append(a)

        divisor = gcd(a, n)
        if divisor != 1:
            res.is_prime = False
            res.reason = f"gcd({a}, {n}) = {divisor} is greater than 1"

            return res.is_prime

        r = pow(a, n - 1, n)

        if r != 1:
            res.is_prime = False
            res.reason = f"r = {a}^({n} - 1) (mod {n}) = {r} is not equal 1"

            return res.is_prime

        return res.is_prime


class SolovayStrassenTest(PrimalityTestABC):
    @property
    def algorithm_name(self) -> str:
        return 'Solovay--Strassen'

    @staticmethod
    def _jacobi_symbol(a: int, n: int) -> int:
        if n < 3 or (n % 2 == 0):
            raise ValueError('n must be an odd integer greater or equal 3')
        if a < 0 or a >= n:
            raise ValueError('a must lie in [0, n)')
        
        g = 1

        while True:
            if a == 0:
                return 0

            if a == 1:
                return g

            a1 = a
            k = 0
            while (a1 & 1) == 0:
                a1 >>= 1
                k += 1

            if ((k & 1) == 1) and ((n & 7) in (3, 5)):
                s = -1
            else:
                s = 1
            
            if a1 == 1:
                return g * s

            if ((n & 3) == 3) and ((a1 & 3) == 3):
                s = -s

            a = n % a1
            n = a1
            g = g * s


    def _perform_round(self, rng: Random, res: TestResult) -> bool:
        n = res.candidate
        a = rng.randint(2, n - 2)

        res.witnesses.append(a)

        symbol = self._jacobi_symbol(a, n)

        r = pow(a, (n - 1) // 2, n)

        if gcd(a, n) != 1 or r != symbol % n:
            res.is_prime = False
            res.reason = f"Euler criterion failed for base {a}: residue {r}, Jacobi symbol {symbol}"

            return res.is_prime

        return res.is_prime


class MillerRabinTest(PrimalityTestABC):
    @property
    def algorithm_name(self) -> str:
        return 'Miller--Rabin'
    
    def _perform_round(self, rng: Random, res: TestResult) -> bool:
        n = res.candidate

        s = 0
        r = n - 1
        while (r & 1) == 0:
            r >>= 1
            s += 1

        a = rng.randint(2, n - 2)
        
        res.witnesses.append(a)
        y = pow(a, r, n)
        if y == 1 or y == n - 1:
            return res.is_prime

        for _ in range(s - 1):
            y = pow(y, 2, n)
            if y == n - 1:
                return res.is_prime
            if y == 1:
                break

        res.is_prime = False
        res.reason = f"Base {a} is a Miller-Rabin witness for compositeness"
        return res.is_prime


@dataclass
class TestResult:
    algorithm: str
    candidate: int
    is_prime: bool = True
    rounds_executed: int = 0
    reason: str = 'Test passed'
    witnesses: deque[int] = field(default_factory=lambda: deque(maxlen=MAX_WITNESSES))


@dataclass
class TestResultPrinter:
    is_verbose: bool

    def _extract_verdict(self, result: TestResult) -> str:
        return 'Possibly prime' if result.is_prime else 'Composite'
    
    def print(self, result: TestResult) -> int:
        if self.is_verbose:
            print('=== RESULT ===')

            print('Algorithm used:', result.algorithm)
            print('n =', result.candidate)

            print('Verdict:', self._extract_verdict(result))

            print('Reason:', result.reason)
            print('Rounds executed:', result.rounds_executed)

            print('Witnesses:', end=' ')
            print(*result.witnesses, sep=', ')

            print()
        else:
            print(
                result.algorithm,
                result.candidate, '->', self._extract_verdict(result)
            )
        
        return 0


class CarmichaelProvider:
    min: int
    mrt: MillerRabinTest

    def __init__(self, min_val: int = 10**29, mrpt: MillerRabinTest | None = None) -> None:
        self.mrt = mrpt if mrpt else MillerRabinTest(15)

        if min_val < 561:
            min_val = 561
        self.min = min_val

    def get_number(self, rng: Random | None = None) -> int:
        if rng is None:
            rng = Random()

        # Calculate minimum k for Chernick's method: (6k+1)(12k+1)(18k+1) >= self.min
        min_k = max(1, isqrt(isqrt(self.min // 1296)) ** 2)
        while (6 * min_k + 1) * (12 * min_k + 1) * (18 * min_k + 1) < self.min:
            min_k += 1

        # Set search window to 10% of min_k (at least 1000) to ensure randomness
        span = max(1000, min_k // 10)

        while True:
            k = rng.randint(min_k, min_k + span)
            p1 = 6 * k + 1
            p2 = 12 * k + 1
            p3 = 18 * k + 1

            # Check if all three factors are prime
            if all([self.mrt.is_prime(p).is_prime for p in {p1, p2, p3}]):
                return p1 * p2 * p3
