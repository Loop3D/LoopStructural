{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}
  
   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module Attributes

   .. autosummary::
      :toctree:                                          
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% set functions_ours = ours(fullname, functions) %}
   {% if functions_ours %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:
      :recursive:
   {% for item in functions_ours %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% set classes_ours = ours(fullname, classes) %}
   {% if classes_ours %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:                                          
      :template: custom-class-template.rst
      :recursive:              
   {% for item in classes_ours %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% set exceptions_ours = ours(fullname, exceptions) %}
   {% if exceptions_ours %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:
      :template: custom-class-template.rst
   {% for item in exceptions_ours %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst                 
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
