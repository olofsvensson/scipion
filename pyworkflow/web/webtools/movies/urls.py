import os
from django.conf.urls import url
import pyworkflow as pw

MEDIA_MOVIES = os.path.join(pw.HOME, 'web', 'webtools', 'movies', 'resources')

urls = [
    (r'^resources_movies/(?P<path>.*)$', 
        'django.views.static.serve', 
        {'document_root': MEDIA_MOVIES}
    ),
    
    url(r'^movies/', 'app.views_webtools.service_movies'),
#     url(r'^mov/', 'app.views_webtools.service_projects'),
    url(r'^check_m_id/$', 'app.views_webtools.check_m_id'),
    url(r'^create_movies_project/$', 'app.views_webtools.create_movies_project'),
#     url(r'^get_testdata/$', 'app.views_webtools.get_testdata'),
    url(r'^m_content/$', 'app.views_webtools.movies_content')
    
]