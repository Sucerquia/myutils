import inspect
from importlib import import_module


def include_doc():
    def decorator(func):
        func.__doc__ += "funcionaaaa"
        return func
    return decorator


# add2executable
def methods_in_class(module, class_name):
    module = import_module(module)

    imported_class = getattr(module, class_name)

    # Get all methods of the class
    methods = [member[0] for member
               in inspect.getmembers(imported_class,
                                     predicate=inspect.isfunction)
               if member[0] != '__init__']
    for method in methods:
        print(method)
    return methods
