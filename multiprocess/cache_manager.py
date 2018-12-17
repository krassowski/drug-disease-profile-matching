from multiprocessing import Manager
from types import FunctionType
from typing import Any, List, Dict

from dataclasses import dataclass


@dataclass
class MultiprocessCache:
    register: Dict[str, FunctionType]
    name: str
    cache_type: str
    instance: Any = None

    def bind_to_register(self):
        register = self.register
        register[self.name] = self.instance

    def create_in(self, manager: Manager):
        constructor = getattr(manager, self.cache_type)
        self.instance = constructor()
        return self.instance


class MultiprocessCacheManager:

    def __init__(self):
        self.caches: List[MultiprocessCache] = []

    def add_cache(self, register, name, cache_type):
        cache = MultiprocessCache(register, name, cache_type)
        self.caches.append(cache)

    def any_cache(self):
        return self.caches[0]

    def respawn_cache_if_needed(self, force=False):
        try:
            for cache in self.caches:
                # TODO: this assumes type=dict
                cache.instance.keys()
                cache.bind_to_register()
            if force:
                raise Exception
        except:
            manager = Manager()
            for cache in self.caches:
                cache.create_in(manager)
                cache.bind_to_register()


multiprocess_cache_manager = MultiprocessCacheManager()
