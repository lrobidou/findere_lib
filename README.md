# findere_lib
header only library of findere


[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

`findere` is a simple strategy for speeding up queries and for reducing false positive calls from any Approximate Membership Query data structure (AMQ). With no drawbacks, queries are two times faster with two orders of magnitudes less false positive calls.

<img src="paper_companion/HMP_analyses/fpr_wrt_z/fpr_simple.png" alt="fpr" width="600"/>

The `findere` library proposed here gives a abstract class, called customAMQ, to wrap your own AMQ and apply findere on it.

For an exemple on how to use findere_lib, see https://github.com/lrobidou/findere.