from abc import abstractmethod, ABCMeta as NativeABCMeta


abstract_method = abstractmethod


def abstract_property(function):
    return property(abstractmethod(function))


class DummyAttribute:
    pass


def abstract_attribute(obj=None):
    if obj is None:
        obj = DummyAttribute()
    obj.__is_abstract_attribute__ = True
    return obj


class ABCMeta(NativeABCMeta):

    def __call__(cls, *args, **kwargs):
        instance = NativeABCMeta.__call__(cls, *args, **kwargs)
        abstract_attributes = {
            name
            for name in dir(instance)
            if getattr(getattr(instance, name), '__is_abstract_attribute__', False)
        }
        if abstract_attributes:
            raise NotImplementedError(
                "Can't instantiate abstract class {} with"
                " abstract attributes: {}".format(
                    cls.__name__,
                    ', '.join(abstract_attributes)
                )
            )
        return instance


class ABC(metaclass=ABCMeta):
    pass


__all__ = ['abstract_method', 'abstract_property', 'abstract_attribute', 'ABCMeta', 'ABC']
