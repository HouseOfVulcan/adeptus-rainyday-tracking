"""User extension point for custom tracker implementations."""


def build_tracker(config):
    """Return a tracker object from user-defined logic."""
    raise NotImplementedError("Implement your custom tracker here")
