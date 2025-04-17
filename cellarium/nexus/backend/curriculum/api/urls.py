from rest_framework.routers import DefaultRouter

from cellarium.nexus.backend.curriculum.api.views import CurriculumViewSet

router = DefaultRouter()
router.register("curriculums", CurriculumViewSet, basename="curriculum")

urlpatterns = router.urls
