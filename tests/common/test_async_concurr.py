'''
Test src/bioomics/common/async_concurr.py
'''
from ..helper import *
import asyncio

from src.bioomics.common.async_concurr import main


async def _double(x):
    await asyncio.sleep(0)
    return x * 2


class TestAsyncConcurr(TestCase):

    def test_main_collects_results(self):
        collected = []
        asyncio.run(main([_double], [1, 2, 3], post_func=collected.append))
        self.assertEqual(collected, [2, 4, 6])

    def test_main_multiple_workers(self):
        collected = []
        asyncio.run(main([_double, _double], [1, 2, 3, 4],
                         post_func=collected.append))
        self.assertEqual(sorted(collected), [2, 4, 6, 8])

    def test_main_without_post_func(self):
        # results are consumed even without a post processor
        asyncio.run(main([_double], [1, 2]))

    def test_main_empty_data(self):
        collected = []
        asyncio.run(main([_double], [], post_func=collected.append))
        self.assertEqual(collected, [])
