from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class BaseBuilder:
    def __init__(self, model, name: str = "Feature"):
        """Base builder that provides a template for
        implementing different builders.


        Parameters
        ----------
        name : str, optional
            The name of the feature being built. The name needs to be unique
            for the model, by default "Feature"

        Notes
        ------
        The interpolation/work should only be done when .build() is called
        .feature should return a reference to the feature, which may not be up
        to date.
        If the build arguments are changed, this will flag that the feature needs to be rebuilt
        """
        self._name = name
        self._model = model
        self._feature = None
        self._up_to_date = False
        self._build_arguments = {}
        self.faults = []

    def set_not_up_to_date(self, caller):
        logger.info(
            f"Setting {self.name} to not up to date from an instance of {caller.__class__.__name__}"
        )
        self._up_to_date = False

    @property
    def model(self):
        return self._model

    @property
    def feature(self):
        return self._feature

    @property
    def build_arguments(self):
        return self._build_arguments

    def update_build_arguments(self, build_arguments):
        # self._build_arguments = {}
        logger.info(f"Setting build arguments for {self.name}")
        for k, i in build_arguments.items():
            logger.info(f"Setting {k} to {i} for {self.name}")
            logger.info(f"{k} is currrently {self._build_arguments.get(k, None)}")
            if i != self._build_arguments.get(k, None):
                logger.info(f"Setting {k} to {i} for {self.name}")
                self._build_arguments[k] = i
                ## if build_arguments change then flag to reinterpolate
                self._up_to_date = False

    def update(self):
        if self._up_to_date:
            logger.info(f"{self.name} is up to date")
            return
        logger.info(f"Updating {self.name}")
        self._feature.faults = self.faults
        self.build(**self.build_arguments)

    def build(self, **kwargs):
        raise NotImplementedError("BaseBuilder should be inherited and build method overwritten")

    @property
    def name(self):
        return self._name

    def up_to_date(self, callback=None):
        """
        check if the feature is uptodate
        if its not update.

        Parameters
        ----------
        callback : function
            a function that is called when the feature is updated

        """
        logger.info(f'Feature: {self.name} up to date: {self._up_to_date}')
        for f in self.faults:
            f.builder.up_to_date(callback=callback)
        # has anything changed in the builder since we built the feature? if so update
        if not self._up_to_date:
            self.update()
            if callable(callback):
                callback(1)
            return
        # check if the interpolator is up to date, if not solve
        if not self._interpolator.up_to_date:
            self.update()
            if callable(callback):
                callback(1)
            return
        if callable(callback):
            callback(1)

    def add_fault(self, fault):
        """
        Add a fault to the geological feature builder

        Parameters
        ----------
        fault : FaultSegment
            A faultsegment to add to the geological feature

        Returns
        -------

        """
        logger.info(f'Adding fault {fault.name} to {self.name}')
        self._up_to_date = False
        self.faults.append(fault)
