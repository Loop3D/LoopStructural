class BaseBuilder:
    def __init__(self, name="Feature"):
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
        self._feature = None
        self._up_to_date = False
        self._build_arguments = {}
        self.faults = []

    @property
    def feature(self):
        return self._feature

    @property
    def build_arguments(self):
        return self._build_arguments

    @build_arguments.setter
    def build_arguments(self, build_arguments):
        # self._build_arguments = {}
        for k, i in build_arguments.items():
            if i != self._build_arguments.get(k, None):
                self._build_arguments[k] = i
                ## if build_arguments change then flag to reinterpolate
                self._up_to_date = False

    def update(self):
        self.build(**self.build_arguments)

    def build(self, **kwargs):
        raise NotImplementedError(
            "BaseBuilder should be inherited and build method overwritten"
        )

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
        for f in self.faults:
            f.builder.up_to_date(callback=callback)
        # has anything changed in the builder since we built the feature? if so update
        if self._up_to_date == False:
            self.update()
            if callable(callback):
                callback(1)
            return
        # check if the interpolator is up to date, if not solve
        if self._interpolator.up_to_date == False:
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
        self._up_to_date = False
        self.faults.append(fault)
