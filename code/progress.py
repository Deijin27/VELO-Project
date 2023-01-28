""" Command line progress bar. """

class ProgressBar:
    """ Command line progress bar controller """

    def __init__(self, total, zero=0, prefix="", suffix="", length=100, decimals=0, fill='█', unfill=" "):
        """ Generate a progress bar and print the first iteration of it.

        Warning - does not work in idle

        Parameters
        ----------
        total : int or float
            The number corresponding to 100% completion.
        zero : int or float, optional
            The number corresponding to 0% completion.
        prefix : str, optional
            The prefix to the progress bar.
        suffix : str, optional
            The suffix to the progress bar.
        length : int, optional
            The number of fill intervals corresponding to 100%.
        decimals : int
            The number of decimals on the percentage progress to be shown.
        fill : str
            The character used to form the filled part of the progress bar.
        unfill : str
            The character used to form the unfilled part of the progres bar.
        """
        self.current = zero

        self.total = total
        self.zero = zero
        self.prefix = prefix
        self.suffix = suffix
        self.length = length
        self.decimals = decimals
        self.fill = fill
        self.unfill = unfill

        self.update(self.current)

    def update(self, current):
        """ Update the progress bar to show an updated percentage.

        Progress
        --------
        current : float or int
            The current progress for the value.
        """
        self.current = current
        percentage = self.current_percentage()
        intervals = int(percentage / 100 * self.length)

        bar = intervals * self.fill + (self.length - intervals) * self.unfill
        str_percentage = format(percentage, f'.{self.decimals}f')

        if not self.current == self.total:
            print_end = '\r'
        else:
            print_end = '\n'

        print(f'\r{self.prefix}|{bar}| {str_percentage}% {self.suffix}', end=print_end)

    def current_percentage(self):
        """ Yield the current percentage of progress. """
        return 100 * (self.current - self.zero) / (self.total - self.zero)
