"""Cached-background drawing with a full-redraw fallback for Matplotlib."""


class BlitRenderer:
    def __init__(self, figure, artists):
        self.figure = figure
        self.canvas = figure.canvas
        self.artist_provider = artists
        self.enabled = bool(self.canvas.supports_blit)
        self.artists = []
        self.background = None
        self.canvas.mpl_connect("draw_event", self._on_draw)
        self.canvas.mpl_connect("resize_event", self.invalidate)
        self.sync()

    def sync(self):
        artists = list(self.artist_provider())
        if artists != self.artists:
            for artist in self.artists:
                artist.set_animated(False)
            self.artists = artists
            for artist in self.artists:
                artist.set_animated(self.enabled)
            self.invalidate()

    def invalidate(self, event=None):
        self.background = None

    def _paint(self):
        for artist in sorted(self.artists, key=lambda a: a.get_zorder()):
            if artist.get_visible() and artist.get_figure() is self.figure:
                self.figure.draw_artist(artist)

    def _on_draw(self, event):
        if not self.enabled or self.canvas.is_saving():
            return
        self.background = self.canvas.copy_from_bbox(self.figure.bbox)
        self._paint()

    def draw(self):
        self.sync()
        if not self.enabled:
            self.canvas.draw_idle()
        elif self.background is None:
            self.canvas.draw()
        else:
            self.canvas.restore_region(self.background)
            self._paint()
            self.canvas.blit(self.figure.bbox)

    def savefig(self, path, **kwargs):
        # Figure-level animated text is otherwise omitted by savefig, even
        # though animated artists inside axes are included when saving.
        self.sync()
        for artist in self.artists:
            artist.set_animated(False)
        try:
            self.figure.savefig(path, **kwargs)
        finally:
            for artist in self.artists:
                artist.set_animated(self.enabled)
            self.invalidate()
