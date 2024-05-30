def include_doc():
    def decorator(func):
        func.__doc__ += "funcionaaaa"
        return func
    return decorator
